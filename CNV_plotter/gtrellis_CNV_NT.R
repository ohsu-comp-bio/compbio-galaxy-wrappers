
VERSION <- "2.0.0"

#### Package Install ####
gtrellisPackageExists <- require ("gtrellis")
if ( !gtrellisPackageExists ) {
	install.packages ("gtrellis")
	library ("gtrellis")
}

fuzzyjoinPackageExists <- require ("fuzzyjoin")
if ( !fuzzyjoinPackageExists ) {
	install.packages("fuzzyjoin")
	library(fuzzyjoin)
}

dplyrPackageExists <- require ("dplyr")
if ( !dplyrPackageExists ) {
	install.packages("dplyr")
	library(dplyr)
}


#########################

#Get the inputs
#send command line arguments to args
args <- commandArgs(TRUE)
#The name of the file containing the read count data
countFile <- args[1]
# the name of the file containing the segment data
segmentFile <- args[2]
#What to name the sample in teh plots
SampleName <- args[3]
#the file containing the relevant genes and thier locations, for coloring the dots
geneFile <- args[4]
#the file containing the chromosome lengths. can be either basic CHR, START, END or a sequence dictionary, for placing the pink lines
chromFile <- args[5]
#get the tumor percent, as a value from 0-100 and convert to fraction
tumorPercent <- as.double(args[6])/100
#the upper and lower bounds of normal to place the pink range lines
upper <- as.double(args[7])
lower <- as.double(args[8])


#find out what type of chrom file was entered by taking the last four characters
chromType = substr(chromFile, nchar(chromFile)-3, nchar(chromFile))

#if the type of file is a sequence dictionary it will end in "dict"
if (chromType == "dict") {
	#if using a sequence dictionary reformat it for use.
	#read in the sequence dictionary, fill in missing values with blanks, there is no header, tab separated
	copiesUpperLower<-read.table(chromFile, fill=TRUE, header = FALSE, sep = "\t")
	#get the chromosome column by splitting the second column by : and taking the second part of the split
	copiesUpperLower$CHR = data.frame(do.call('rbind', strsplit(as.character(copiesUpperLower[,2]), ":", fixed=TRUE)))$X2
	#All the starts are 1, make teh start column
	copiesUpperLower$START = 1
	#Get the end cloumn usign the same method as teh chr column but splitting the third column instead.
	copiesUpperLower$END = data.frame(do.call('rbind', strsplit(as.character(copiesUpperLower[,3]), ":", fixed=TRUE)))$X2
	#remove all irrelevant columns
	copiesUpperLower = copiesUpperLower[,c("CHR", "START", "END")]
	#properly format the columns as character or numbers
	copiesUpperLower$CHR <- lapply(copiesUpperLower$CHR, as.character)
	copiesUpperLower$END <- lapply(copiesUpperLower$END, as.character)
	copiesUpperLower$END <- sapply(copiesUpperLower$END, as.numeric)
	#remove all rows that have non-chromosome contigs
	copiesUpperLower = subset(copiesUpperLower, nchar(CHR) < 3)
	#remove the mitochondrial contig
	copiesUpperLower = subset(copiesUpperLower, CHR != "MT")
	#add the required "chr" to the contig names
	copiesUpperLower$CHR = paste0("chr", copiesUpperLower$CHR)
} else {
	#if using the already formatted chrom file simply read it in wiht a header. should have CHR, START and END, headers and be separated by tabs
	copiesUpperLower<-read.delim(chromFile, header = TRUE, sep = "\t")
}


#read in the files for the read counts, segemnts and genefiel.
#read counts will create the blue dots
countsData <- read.table(countFile, header = TRUE, sep = "\t", comment.char = "@")
#segements create teh green lines
segmentsData <- read.table(segmentFile, header = TRUE, sep = "\t", comment.char = "@")
#gene intervals color the dots based on revelance
geneIntervals <- read.table(geneFile, header = TRUE, sep = "\t", comment.char = "@")

#a function to get the correct number of copies based on teh tumor percent. Pass it the data sheet in question and which column has the uncorrected copy number
correctCopies <- function(data_sheet) {
	#remove empty lines
	data_sheet<-data.frame(data_sheet[complete.cases(data_sheet),])
	log_name <- colnames(data_sheet[endsWith(colnames(data_sheet), "LOG2_COPY_RATIO")])
	#get the tumor corrected copies number
	data_sheet$tumorCorrectedCopies = (((2^(data_sheet[,log_name] + 1))-(2*(1-tumorPercent)))/tumorPercent)
	#replace the negatives with .1
	data_sheet$Tumor_Corrected_Copies_STPv3 = (ifelse(data_sheet$tumorCorrectedCopies > 0.1, data_sheet$tumorCorrectedCopies, .1))
	#get the raw copy number
	data_sheet$rawCopies = 2*(2^data_sheet[,log_name])
	#return only the needed columns, contig start end and tumore corrected copies
	return(data_sheet[,c("CONTIG", "START", "END", log_name, "rawCopies", "Tumor_Corrected_Copies_STPv3")])
}


#Get the tumor corrected copy number for both the coutns and the segments
cnvPoints <- correctCopies(countsData)
cnvSegments <- correctCopies(segmentsData)


#get the stpv3_225 data from the gene intervals data
#create a new column with the cromosomal location in a single column for mergeing with the gene interval file
cnvPoints$hg19Loc = paste0("chr", cnvPoints$CONTIG, ":", cnvPoints$START, "-", cnvPoints$END)
#merge the gene interval information with the cnv points information
cnvPoints <- merge(cnvPoints, geneIntervals)
#remove empty lines
cnvPoints<-data.frame(cnvPoints[complete.cases(cnvPoints),])

#add the needed "chr" to each chromosome name in the counts dataset and teh segments dataset
cnvPoints$CONTIG = paste0("chr", cnvPoints$CONTIG)
cnvSegments$CONTIG = paste0("chr", cnvSegments$CONTIG)

## Tumor-corrected counts per segment: change the tumor content adjustment of the individual counts per segment.
cnvPointsSegments <- genome_full_join(cnvPoints, cnvSegments, c('CONTIG', 'START', 'END'))
cnvPointsSegments<-data.frame(cnvPointsSegments[complete.cases(cnvPointsSegments),]) #? merge() introducing NAs in some samples
colnames(cnvPointsSegments) <- gsub(".x$", "_point", colnames(cnvPointsSegments))
colnames(cnvPointsSegments) <- gsub(".y$", "_segment", colnames(cnvPointsSegments))
cnvPointsSegments$tumorCorrectedCopiesPerSegment <- cnvPointsSegments$rawCopies_point * (cnvPointsSegments$Tumor_Corrected_Copies_STPv3_segment / cnvPointsSegments$rawCopies_segment)
cnvPointsSegments$hg19Loc_segment = paste0(cnvPointsSegments$CONTIG_segment, ":", cnvPointsSegments$START_segment, "-", cnvPointsSegments$END_segment)
segment_log2rawCopies_mad <- cnvPointsSegments %>%
  group_by(hg19Loc_segment) %>%
  mutate(log2rawCopies_mad_segment = mad(LOG2_COPY_RATIO)) %>%
  select(hg19Loc_segment,
         log2rawCopies_mad_segment) %>%
  unique()

sample_log2rawCopies_mad <- round(median(segment_log2rawCopies_mad$log2rawCopies_mad_segment), digits = 2)
print(paste0('Sample Median MAD ("Median of the segment-level MADD values of log2 raw copies per sample"): ', sample_log2rawCopies_mad))
write(paste0('{{\"cnv_median_segment_mad_cn\": ', sample_log2rawCopies_mad, '}}'), "output_metrics.txt")


plotCNV <- function(data, copyratioCol, copyratioPerSegmentCol, panelGeneCol, segments,
					filename, SampleName, title, ncol, byrow, xaxis, add_ideogram_track) {

	png(filename = paste0(filename, ".png"), width = 1920, height = 960, units = "px")
	gtrellis_layout(n_track = 1, title = paste0(SampleName, ': ', title), title_fontsize = 28, ncol = ncol, byrow = byrow,
                	track_axis = TRUE, track_ylim = c(-5, (log2(max(data[,paste0(copyratioCol, '_point')])) + 1)),
                	lab_fontsize = 20, axis_label_fontsize = 16, xaxis = xaxis, add_name_track = TRUE, name_fontsize = 16, add_ideogram_track = add_ideogram_track)
	# color points differently if targets are within or outside STPv3 genes of interest
	# blue == within STPv3 gene of interest; grey == outside STPv3 gene of interest; red == neither
	add_points_track(data, log2(data[,copyratioPerSegmentCol]), pch = 16, size = unit(2, "mm"),
					 gp = gpar(col = ifelse(data[,panelGeneCol] == "1", "blue",
                                 ifelse(data[,panelGeneCol] == "0", "grey", "red"))))
	#add the green segement lines
	add_segments_track(segments, log2(segments[,copyratioCol]), track=current_track(), gp = gpar(col = "green", lwd = 5))

	# add lines across plots to demarcate 5 and 0.5 tumor-corrected copies; log2(5)=2.32; log2(0.5)= -1
	# the location of the lines at 5 and .5 are now moveable with the upper and lower arguments.
	add_segments_track(copiesUpperLower, log2(upper),track=current_track(), gp = gpar(col = "pink", lty = 2, lwd = 1))
	add_segments_track(copiesUpperLower, log2(lower),track=current_track(), gp = gpar(col = "pink", lty = 2, lwd = 1))

	# complete png export plot1
	dev.off()
}

## First plot of all chr in one row
plotCNV(cnvPointsSegments, "Tumor_Corrected_Copies_STPv3", "tumorCorrectedCopiesPerSegment", "STPv3_225", cnvSegments,
		"plot1", SampleName, "Log2 tumor-corrected copies", 24, TRUE, FALSE, FALSE)
## Second plot of all chr in six columns
# plot2 has six columns; chr ordered by column, not by row
# ylim max is set to max segment value for sample, plus add 1 for buffer to prevent points from getting squished at top
# ylim min is fixed at -5 for all samples; may cut out failed assays, but prevents individual failed assays from making min too low
plotCNV(cnvPointsSegments, "Tumor_Corrected_Copies_STPv3", "tumorCorrectedCopiesPerSegment", "STPv3_225", cnvSegments,
		"plot2", SampleName, "Log2 tumor-corrected copies", 6, FALSE, TRUE, TRUE)

#### Duplicate plot1 and plot2 to make additional plots with raw copies ####
#### relabel plot1 to plot3 and relabel plot2 to plot4 for the raw copies plots ####
#### change plot titles, track_ylim, and track_ylab from "tumor-corrected copies" to "raw copies" ####
#### Replace cnvPoints$Tumor_Corrected_Copies_STPv3 with cnvPoints$rawCopies #####
#### Replace cnvSegments$Tumor_Corrected_Copies_STPv3 with cnvSegments$rawCopies #####
plotCNV(cnvPointsSegments, "rawCopies", "tumorCorrectedCopiesPerSegment", "STPv3_225", cnvSegments,
		"plot3", SampleName, "Log2 raw copies", 24, TRUE, FALSE, FALSE)
plotCNV(cnvPointsSegments, "rawCopies", "tumorCorrectedCopiesPerSegment", "STPv3_225", cnvSegments,
		"plot4", SampleName, "Log2 raw copies", 6, FALSE, TRUE, TRUE)


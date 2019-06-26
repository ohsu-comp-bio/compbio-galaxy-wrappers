

#### Package Install ####
gtrellisPackageExists <- require ("gtrellis")

if ( !gtrellisPackageExists ) {
  install.packages ("gtrellis")
  library ("gtrellis")
  
}

#########################

#Get the inputs
#send command line arguments to args
args <- commandArgs(TRUE)
#get the tumor percent, as a value from 0-100
tumorPercent <- as.double(args[1])
#The name of the file containing the read count data
countFile <- args[2]
# the name of the file containing the segment data
segmentFile <- args[3]
#What to name the sample in teh plots
SampleName <- args[4]
#the file containing the relevant genes and thier locations, for coloring the dots
geneFile <- args[5]
#the file containing the chromosome lengths. can be either basic CHR, START, END or a sequence dictionary, for placing the pink lines
chromFile <- args[6]
#the upper and lower bounds of normal to place the pink range lines
upper <- as.double(args[7])
lower <- as.double(args[8])

#turn tumor percent into a fraction
tumorPercent = tumorPercent/100

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
correctCopies <- function(data_sheet, log_name) {
	#remove empty lines
	data_sheet<-data.frame(data_sheet[complete.cases(data_sheet),])
	#get the tumor corrected copies number
	data_sheet$tumorCorrectedCopies = (((2^(log_name + 1))-(2*(1-tumorPercent)))/tumorPercent)
	#replace the negatives with .1
	data_sheet$Tumor_Corrected_Copies_STPv3 = (ifelse(data_sheet$tumorCorrectedCopies > 0.1, data_sheet$tumorCorrectedCopies, .1))
	#return only the needed columns, contig start end and tumore corrected copies
	return(data_sheet[,c("CONTIG", "START", "END", "Tumor_Corrected_Copies_STPv3")])
}

#Get the tumor corrected copy number for both the coutns and the segments
cnvPoints <- correctCopies(countsData, countsData$LOG2_COPY_RATIO)
cnvSegments <- correctCopies(segmentsData, segmentsData$MEAN_LOG2_COPY_RATIO)

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

## First plot of all chr in one row  
# create png for plot1
png(filename = "plot1.png", 
    width = 1920, 
    height = 960, 
    units = "px")


# plot1 has one row for all chr 
# ylim max is set to max point value for sample, plus add 1 for buffer to prevent points from getting squished at top
# ylim min is fixed at -5 for all samples; may cut out failed assays, but prevents individual failed assays from making min too low
gtrellis_layout(n_track = 1, 
                title = SampleName,
                title_fontsize = 32,
                track_axis = TRUE,
                track_ylim = c(-5, (log2(max(cnvPoints$Tumor_Corrected_Copies_STPv3)) + 1)),
                track_ylab = "Log2 tumor-corrected copies",
                lab_fontsize = 20,
                axis_label_fontsize = 18,
                xaxis = FALSE,
                add_name_track = TRUE,
                name_fontsize = 14,
                add_ideogram_track = FALSE) 



# color points differently if targets are within or outside STPv3 genes of interest
# blue == within STPv3 gene of interest; grey == outside STPv3 gene of interest; red == neither

add_points_track(cnvPoints, log2(cnvPoints$Tumor_Corrected_Copies_STPv3),
                 pch = 16, 
                 size = unit(2, "mm"),
                 gp = gpar(col = ifelse(cnvPoints$STPv3_225 == "1", "blue",
                                 ifelse(cnvPoints$STPv3_225 == "0", "grey", "red")))
                ) 

#add the green segement lines
add_segments_track(cnvSegments, log2(cnvSegments$Tumor_Corrected_Copies_STPv3),track=current_track(),
                   gp = gpar(col = "green",
                             lwd = 5)
                   ) 
   
# add lines across plots to demarcate 5 and 0.5 tumor-corrected copies; log2(5)=2.32; log2(0.5)= -1
# the location of the lines at 5 and .5 are now moveable with the upper and lower arguments.              
add_segments_track(copiesUpperLower, log2(upper),track=current_track(),
                gp = gpar(col = "pink",
                lty = 2,
                lwd = 1)
                )
                
add_segments_track(copiesUpperLower, log2(lower),track=current_track(),
                gp = gpar(col = "pink",
                lty = 2,
                lwd = 1)
                )                   

# complete png export plot1
dev.off()


  
## Second plot of all chr in six columns 
# create png for plot2
png(filename = "plot2.png",
    width = 1920,
    height = 960,
    units = "px")



# plot2 has six columns; chr ordered by column, not by row
# ylim max is set to max segment value for sample, plus add 1 for buffer to prevent points from getting squished at top
# ylim min is fixed at -5 for all samples; may cut out failed assays, but prevents individual failed assays from making min too low
gtrellis_layout(n_track = 1, 
                title = paste0(SampleName, ": Log2 tumor-corrected copies",sep=""),
                title_fontsize = 32,
                ncol = 6,
                byrow = FALSE,
                track_axis = TRUE,
                track_ylim = c(-5, (log2(max(cnvPoints$Tumor_Corrected_Copies_STPv3)) + 1)),
                lab_fontsize = 20,
                axis_label_fontsize = 16,
                add_name_track = TRUE,
                name_fontsize = 18,
                add_ideogram_track = TRUE) 


# color points differently if targets are within or outside STPv3 genes of interest
# blue == within STPv3 gene of interest; grey == outside STPv3 gene of interest; red == neither

add_points_track(cnvPoints, log2(cnvPoints$Tumor_Corrected_Copies_STPv3),
                 pch = 16, 
                 size = unit(1.5, "mm"),
                 gp = gpar(col = ifelse(cnvPoints$STPv3_225 == "1", "blue",
                                        ifelse(cnvPoints$STPv3_225 == "0", "grey", "red")))
                )  

add_segments_track(cnvSegments, log2(cnvSegments$Tumor_Corrected_Copies_STPv3),track=current_track(),
                   gp = gpar(col = "green",
                             lwd = 2)
                   ) 
# add lines across plots to demarcate 5 and 0.5 tumor-corrected copies; log2(5)=2.32; log2(0.5)= -1
# the location of the lines at 5 and .5 are now moveable with the upper and lower arguments.            
add_segments_track(copiesUpperLower, log2(upper),track=current_track(),
                gp = gpar(col = "pink",
                lty = 2,
                lwd = 1)
                )
                
add_segments_track(copiesUpperLower, log2(lower),track=current_track(),
                gp = gpar(col = "pink",
                lty = 2,
                lwd = 1)
                )

# complete png export plot2
dev.off()

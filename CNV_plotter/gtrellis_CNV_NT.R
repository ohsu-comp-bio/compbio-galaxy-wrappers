

#### Package Install ####
gtrellisPackageExists <- require ("gtrellis")

if ( !gtrellisPackageExists ) {
  install.packages ("gtrellis")
  library ("gtrellis")
  
}

# tidyversePackageExists <- require ("tidyverse")

# if ( !tidyversePackageExists ) {
  # install.packages ("tidyverse")
  # library ("tidyverse")
 # }
#########################

#Get the inputs
args <- commandArgs(TRUE)
tumorPercent <- as.double(args[1])
countFile <- args[2]
segmentFile <- args[3]
SampleName <- args[4]
geneFile <- args[5]
chromFile <- args[6]
upper <- as.double(args[7])
lower <- as.double(args[8])

#turn tumor percent into a fraction

tumorPercent = tumorPercent/100

# read file that contains hg19 chromosome BED, along with max and min copy values values that will be plotted as lines across each plot
# working directory is "C:/Users/beadling/Documents"
copiesUpperLower<-read.delim(chromFile, header = TRUE, sep = "\t")


#read in the files from galaxy steps 100, 111, and 115
countsData <- read.table(countFile, header = TRUE, sep = "\t", comment.char = "@")
segmentsData <- read.table(segmentFile, header = TRUE, sep = "\t", comment.char = "@")
geneIntervals <- read.table(geneFile, header = TRUE, sep = "\t", comment.char = "@")

#geneIntervals <- mutate(geneIntervals, hg19Loc = substring(STRINGhg19, 4))
#geneIntervals$hg19Loc = substr(geneIntervals$STRINGhg19, 4)

correctCopies <- function(data_sheet, log_name) {
	#remove empty lines
	data_sheet<-data.frame(data_sheet[complete.cases(data_sheet),])
	#get the tumor corrected copies number
	#data_with_neg <- mutate(data_sheet, tumorCorrectedCopies = ((2^(log_name + 1))-(2*(1-tumorPercent))/tumorPercent))
	data_sheet$tumorCorrectedCopies = ((2^(log_name + 1))-(2*(1-tumorPercent))/tumorPercent)
	#replace the negatives with .1
	#data_no_neg <- mutate(data_with_neg, Tumor_Corrected_Copies_STPv3 = if_else(tumorCorrectedCopies > 0, tumorCorrectedCopies, 0.1))
	data_sheet$Tumor_Corrected_Copies_STPv3 = (ifelse(data_sheet$tumorCorrectedCopies > 0, data_sheet$tumorCorrectedCopies, .1))
	#data_sheet$tumorCorrectedCopies[data_sheet$tumorCorrectedCopies < 0] <- .1
	#data_sheet$Tumor_Corrected_Copies_STPv3 = data_sheet$tumorCorrectedCopies
	#return only the needed columns, contig to end and tumore corrected copies
	#return(select(data_no_neg, CONTIG:END, Tumor_Corrected_Copies_STPv3))
	return(data_sheet[,c("CONTIG", "START", "END", "Tumor_Corrected_Copies_STPv3")])
}

#do the mutations on both data sets
cnvPoints <- correctCopies(countsData, countsData$LOG2_COPY_RATIO)
cnvSegments <- correctCopies(segmentsData, segmentsData$MEAN_LOG2_COPY_RATIO)

#get the stpv3_225 data from the gene intervals data
#create a single location column
#cnvPoints <- mutate(cnvPoints, locationStart = paste(CONTIG,START, sep=":"))
#cnvPoints <- mutate(cnvPoints, hg19Loc = paste(locationStart,END, sep="-"))
cnvPoints$hg19Loc = paste0("chr", cnvPoints$CONTIG, ":", cnvPoints$START, "-", cnvPoints$END)
#merge the gene interval information with the cnv points information
#print.table(cnvPoints)
#print.table(geneIntervals)
cnvPoints <- merge(cnvPoints, geneIntervals)
#remove empty lines
#print.table(cnvPoints)
cnvPoints<-data.frame(cnvPoints[complete.cases(cnvPoints),])

#cnvPoints <- select(cnvPoints, CONTIG, START, END, STPv3_225, Tumor_Corrected_Copies_STPv3)
cnvPoints = cnvPoints[,c("CONTIG", "START", "END", "STPv3_225", "Tumor_Corrected_Copies_STPv3")]
#cnvSegments <- select(cnvSegments, CONTIG, START, END, Tumor_Corrected_Copies_STPv3)
cnvSegments = cnvSegments[,c("CONTIG", "START", "END", "Tumor_Corrected_Copies_STPv3")]

#cnvPoints <- mutate(cnvPoints, CONTIG = paste("chr", CONTIG, sep=''))
cnvPoints$CONTIG = paste0("chr", cnvPoints$CONTIG)
#cnvSegments <- mutate(cnvSegments, CONTIG = paste("chr", CONTIG, sep=''))
cnvSegments$CONTIG = paste0("chr", cnvSegments$CONTIG)	


# # Make list of Excel files that contain CNV data; each file has multiple sheets
# file.list <- list.files(path = "C:/Users/beadling/Desktop/STPv3CNV/", pattern='*_STPv3_CNV.xlsx', full.names = TRUE)

# # Read Excel sheets with CNV segment data; append filename in Source column; combine data for all samples into one file
# df.list.segments <- lapply(file.list, read_excel, sheet = "SEGMENTS")
  # names(df.list.segments)<-file.list
  # for (i in file.list)
  # df.list.segments[[i]]$Source = i
  # do.call(rbind, df.list.segments)  
   
# # Read Excel sheets with CNV counts data; append filename in Source column;  combine data for all samples into one file    
# df.list.counts <- lapply(file.list, read_excel, she~/Box Sync/CNV/gene_intervals.txtet = "COUNTS")
  # names(df.list.counts)<-file.list
  # for (i in file.list)
  # df.list.counts[[i]]$Source = i
  # do.call(rbind, df.list.counts)   
  
# # Create SampleName, segments, and points from each sample Excel file in file.list
# for(i in 1:length(file.list)){
 # SampleName <- basename(file.list[[i]])
 # cnvSegments<-df.list.segments[[file.list[i]]]
 # cnvPoints<-df.list.counts[[file.list[i]]]
 
 

# Low % tumor cases may have negative tumor-corrected copy number, which generates error upon log2 transformation
# Corrected these in Excel prior to import into R, by recoding negative tumor-corrected copy number values to positive number 0.1
# Exclude rows that have missing values
#cnvPoints<-data.frame(cnvPoints[complete.cases(cnvPoints),])
#cnvSegments<-data.frame(cnvSegments[complete.cases(cnvSegments),])
 
 
## First plot of all chr in one row  
# create png for plot1
png(filename = "plot1.png", 
    width = 1920, 
    height = 960, 
    units = "px")


# plot1 has one row for all chr 
# remove .xlsx from end of SampleName
# ylim max is set to max point value for sample, plus add 1 for buffer to prevent points from getting squished at top
# ylim min is fixed at -5 for all samples; may cut out failed assays, but prevents individual failed assays from making min too low
#print.table(cnvPoints)
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
# add lines across plots to demarcate 5 and 0.5 tumor-corrected copies; log2(5)=2.32; log2(0.5)= -1

add_points_track(cnvPoints, log2(cnvPoints$Tumor_Corrected_Copies_STPv3),
                 pch = 16, 
                 size = unit(2, "mm"),
                 gp = gpar(col = ifelse(cnvPoints$STPv3_225 == "1", "blue",
                                 ifelse(cnvPoints$STPv3_225 == "0", "grey", "red")))
                ) 


add_segments_track(cnvSegments, log2(cnvSegments$Tumor_Corrected_Copies_STPv3),track=current_track(),
                   gp = gpar(col = "green",
                             lwd = 5)
                   ) 
                   
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
# add_segments_track(copiesMaxMin, log2(copiesUpperLower$UPPER),track=current_track(),
                # gp = gpar(col = "pink",
                # lty = 2,
                # lwd = 1)
                # )

# add_segments_track(copiesMaxMin, log2(copiesUpperLower$LOWER),track=current_track(),
                   # gp = gpar(col = "pink",
                             # lty = 2,
                             # lwd = 1)
                  # )

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
# add lines across plots to demarcate 5 and 0.5 tumor-corrected copies; log2(5)=2.32; log2(0.5)= -1
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

# add_segments_track(copiesMaxMin, log2(copiesUpperLower$UPPER),track=current_track(),
                   # gp = gpar(col = "pink",
                             # lty = 2,
                             # lwd = 0.1)
                  # )

# add_segments_track(copiesMaxMin, log2(copiesUpperLower$LOWER),track=current_track(),
                   # gp = gpar(col = "pink",
                             # lty = 2,
                             # lwd = 0.1)
                             # )

# complete png export plot2
dev.off()

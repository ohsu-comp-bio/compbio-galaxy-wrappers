

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
countFile <- args[1]
SampleName <- args[2]
chromFile <- args[3]
upper <- as.double(args[4])
lower <- as.double(args[5])



# read file that contains hg19 chromosome BED, along with max and min copy values values that will be plotted as lines across each plot
copiesUpperLower<-read.delim(chromFile, header = TRUE, sep = "\t")

#read in the files
cnvPoints <- read.table(countFile, header = TRUE, sep = "\t", comment.char = "@")

#create a single location column
cnvPoints$hg19Loc = paste0("chr", cnvPoints$Chr, ":", cnvPoints$Start, "-", cnvPoints$End)
cnvPoints = cnvPoints[,c("Chr", "Start", "End", "Copy.Number", "GENE")]
colnames(cnvPoints) <- c("Chr", "Start", "End", "CopyNumber", "GENE")
cnvPoints <- cnvPoints[complete.cases(cnvPoints[,"CopyNumber"]),]

#cnvPoints <- mutate(cnvPoints, CONTIG = paste("chr", CONTIG, sep=''))
cnvPoints$Chr = paste0("chr", cnvPoints$Chr)	

 
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
                track_ylim = c(-1, max(cnvPoints$CopyNumber)),
                track_ylab = "Copy.Number",
                lab_fontsize = 20,
                axis_label_fontsize = 18,
                xaxis = FALSE,
                add_name_track = TRUE,
                name_fontsize = 14,
                add_ideogram_track = FALSE) 



# color points differently if targets are within or outside STPv3 genes of interest
# blue == within STPv3 gene of interest; grey == outside STPv3 gene of interest; red == neither
# add lines across plots to demarcate 5 and 0.5 tumor-corrected copies; log2(5)=2.32; log2(0.5)= -1
add_segments_track(cnvPoints, (cnvPoints$CopyNumber),
                 #pch = 16, 
                 #size = unit(2, "mm"),
                 gp = gpar(col = ifelse(is.na(cnvPoints$GENE), "grey", "blue"))
                ) 
                
add_points_track(cnvPoints, (cnvPoints$CopyNumber),track=current_track(),
                 pch = 16, 
                 size = unit(1.5, "mm"),
                 gp = gpar(col = ifelse(is.na(cnvPoints$GENE), "grey", "blue"))
                 )


add_segments_track(copiesUpperLower, (upper),track=current_track(),
                gp = gpar(col = "pink",
                lty = 2,
                lwd = 1)
                )
                
add_segments_track(copiesUpperLower, (lower),track=current_track(),
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
                title = paste0(SampleName, ": Copy.Number",sep=""),
                title_fontsize = 32,
                ncol = 6,
                byrow = FALSE,
                track_axis = TRUE,
                track_ylim = c(-1, max(cnvPoints$CopyNumber)),
                lab_fontsize = 20,
                axis_label_fontsize = 16,
                add_name_track = TRUE,
                name_fontsize = 18,
                add_ideogram_track = TRUE) 


# color points differently if targets are within or outside STPv3 genes of interest
# blue == within STPv3 gene of interest; grey == outside STPv3 gene of interest; red == neither
# add lines across plots to demarcate 5 and 0.5 tumor-corrected copies; log2(5)=2.32; log2(0.5)= -1

add_segments_track(cnvPoints, (cnvPoints$CopyNumber),
                 #pch = 16, 
                 #size = unit(1.5, "mm"),
                 gp = gpar(col = ifelse(is.na(cnvPoints$GENE), "grey", "blue"))
                 )
                 
add_points_track(cnvPoints, (cnvPoints$CopyNumber),track=current_track(),
                 pch = 16, 
                 size = unit(1.5, "mm"),
                 gp = gpar(col = ifelse(is.na(cnvPoints$GENE), "grey", "blue"))
                 )   
                                 
add_segments_track(copiesUpperLower, (upper),track=current_track(),
                gp = gpar(col = "pink",
                lty = 2,
                lwd = 1)
                )
                
add_segments_track(copiesUpperLower, (lower),track=current_track(),
                gp = gpar(col = "pink",
                lty = 2,
                lwd = 1)
                )
                
 


# complete png export plot2
dev.off()


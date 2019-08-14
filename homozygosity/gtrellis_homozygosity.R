

#### Package Install ####
gtrellisPackageExists <- require ("gtrellis")

if ( !gtrellisPackageExists ) {
  install.packages ("gtrellis")
  library ("gtrellis")
  
}

#########################

#Get the inputs
args <- commandArgs(TRUE)
zygosityFile <- args[1]



zygosityData <- read.table(text = gsub("/", "\t", readLines(zygosityFile)), header = FALSE, skip = 1)


zygosityData$zygosity = (zygosityData$V3==zygosityData$V4)
zygosityData$Homozygous = as.integer(as.logical(zygosityData$zygosity))


colnames(zygosityData) <- c("CHROM", "POS", "GT1", "GT2", "zygosity", "Homozygous")
zygosityData$END = zygosityData$POS + 1
#head(zygosityData, n=10)

# for (chrom in 1:25) {
	# if (chrom == 23) {
		# chrom = "X"
	# } else if (chrom == 24) {
		# chrom = "Y"
	# } else if (chrom == 25) {
		# chrom = "MT"
	# }

# plotData <- subset(zygosityData, CHROM == chrom, select=c(POS, Homozygous))

# png_name = paste0(chrom, "plot.png")

# png(png_name)

# plot(plotData$POS, plotData$Homozygous)

# dev.off()

# }

gData <- subset(zygosityData, select=c(CHROM, POS, END, Homozygous))
colnames(gData) <- c("CHR", "START", "END", "Homozygous")

gData$CHR = paste0("chr", gData$CHR)

#head(gData, n=10)

 
# ## First plot of all chr in one row  
# # create png for plot1
 # png(filename = "plot1.png", 
     # width = 1920, 
     # height = 960, 
     # units = "px")


# # plot1 has one row for all chr 
# gtrellis_layout(n_track = 1, 
                # title = "SampleName",
                # title_fontsize = 32,
                # track_axis = TRUE,
                # track_ylim = c(-1, 2),
                # track_ylab = "Zygosity",
                # lab_fontsize = 20,
                # axis_label_fontsize = 18,
                # xaxis = FALSE,
                # add_name_track = TRUE,
                # name_fontsize = 14,
                # add_ideogram_track = FALSE) 

# add_points_track(gData, gData$Homozygous,
                 # pch = 16, 
                 # size = unit(2, "mm"),
                 # gp = gpar(col = "green", lwd = 5)
                # ) 


# # complete png export plot1
# dev.off()

  
## Second plot of all chr in six columns 
# create png for plot2
png(filename = "plot.png",
    width = 1920,
    height = 960,
    units = "px")



# plot2 has six columns; chr ordered by column, not by row
gtrellis_layout(n_track = 1, 
                title = "Zygosity",
                title_fontsize = 32,
                ncol = 6,
                byrow = FALSE,
                track_axis = TRUE,
                track_ylim = c(-1, 2),
                lab_fontsize = 20,
                axis_label_fontsize = 16,
                add_name_track = TRUE,
                name_fontsize = 18,
                add_ideogram_track = TRUE) 


add_points_track(gData, gData$Homozygous,
                 pch = 16, 
                 size = unit(1.5, "mm"),
                 gp = gpar("blue")  
                 )

# complete png export plot2
dev.off()

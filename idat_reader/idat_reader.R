# VERSION 1.0.1

library(illuminaio)
library(reshape2)
library(tibble)

### START ARGS
args <- commandArgs(trailingOnly=TRUE)
red <- args[1]
reen <- args[2]
manifest <- args[3]
output_filename <- args[4]

# Example
red <-'./idat/Galaxy2-[206112280126_R08C01_Red.idat].idat'
green<-'./idat/Galaxy1-[206112280126_R08C01_Grn.idat].idat'
manifest<-'infinium-methylationepic-v-1-0-b5-manifest-file.csv'

# Read in .idat files
idat_red <- readIDAT(red)
idat_green <- readIDAT(green)

# Organize and merge green and red array data frames
red.counts <- tibble::rownames_to_column(data.frame(idat_red$Quants), "IlmnID")
colnames(red.counts) <- c('IlmnID','red.mean','red.sd','red.NBead')
green.counts <- tibble::rownames_to_column(data.frame(idat_green$Quants), "IlmnID")
colnames(green.counts) <- c('IlmnID','green.mean','green.sd','green.NBead')
total.counts <- merge(red.counts,green.counts,by='IlmnID')
total.counts$IlmnID = paste0('cg', total.counts$IlmnID)

# Load in Illumina Methylation manifest + write to text
mani <- read.csv(manifest,skip=7)
output <- merge(total.counts, mani, by='IlmnID')
write.table(output,output_filename,sep=',',row.names=FALSE)

# VERSION 1.0.0

library(illuminaio)
library(reshape2)
library(tibble)

### START ARGS
args <- commandArgs(trailingOnly=TRUE)
red <- args[1]
green <- args[2]
manifest <- args[3]

# Example
#green<-'Galaxy1-[206137490091_R05C01_Grn.idat].idat'
#manifest<-'infinium-methylationepic-v-1-0-b5-manifest-file.csv'

# Revisit -- use RegEx to use extract specimen ID to name output file
filename<-'idat_output.txt'

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
write.table(output,filename,sep=',',row.names=FALSE)

# Title     : Methylation QC
# Objective : Collect QC metrics that the lab wants to see.
# Created by: letaw
# Created on: 4/2/21\
# Version: 1.0.0

load("mnp.RData")
library(jsonlite)
library(mnp.v11b4)
library(mnpqc)

args <- commandArgs(TRUE)
path <- args[1]
sample <- args[2]
sample_path <- paste(path,sample,sep='/')

# read data using standard minfi method
meth <- read.metharray(sample_path, verbose=FALSE)
# QC plot, mildly useful.
# Output file names look like: MNPqc_blah.png, where blah is the sample name.
mnpqc::MNPqc(meth, detailplot=TRUE)
det <- detectionP(meth)

# Set failed probes.
failed <- det>0.01
# Find p-values for predictive probes.
cg12434587_p <- det['cg12434587',]
cg12981137_p <- det['cg12981137',]
# Total number of probes.
probe_count <- nrow(failed)
# Total number of failed probes.
probe_failed <- sum(rowMeans(failed))
df <- data.frame(probe_count, probe_failed, cg12434587_p, cg12981137_p)

to_write <- jsonlite::toJSON(df[1,], pretty = FALSE, auto_unbox = TRUE)
write(to_write, file = "qc_metrics.json")

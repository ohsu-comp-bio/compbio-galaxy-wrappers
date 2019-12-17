#!/usr/bin/env Rscript

ruvSeqPackageExists <- require ("RUVSeq")
if ( !ruvSeqPackageExists ) {
  install.packages ("RUVSeq")}

library("RUVSeq")

#Get the inputs
args <- commandArgs(TRUE)
dataFile <- args[1]
factors <- as.double(args[2])
minNumReads <- as.double(args[3])
minNumSamples <- as.double(args[4])
repID <- args[5]
cIdx <- args[6]

#read datafile
data <- read.table(dataFile, header = TRUE, sep="\t", row.names=1)

#filter
filter <-  apply(data, 1, function(x) length(x[x>minNumReads])>=minNumSamples)
filtered <- data[filter,]


#genes <- rownames(filtered)[grep(cIdx, rownames(filtered))]
#spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

#select genes to use for normalization
if (cIdx=="all"){
	genes <- rownames(filtered)
} else if (cIdx=="regex"){
        genes = rownames(filtered)[grep(cIdx, rownames(filtered))]
} else if (cIdx=="genes"){
	genes <- strsplit(cIdx, ",")
	genes <- unlist(genes, use.names=FALSE)
} else {
      print("invalid cIdx")
}

#get index of sample replicates
samp_reps <- colnames(filtered)

samp_reps[(grepl(repID, samp_reps))]="replicate_control"

#prepare for ruvseq function
matx <- as.matrix(filtered)
matfiltered <- round(matx, 0)


factorX <- as.factor(c(samp_reps))
groupings <- makeGroups(factorX)

ruv <- RUVs(matfiltered, cIdx=genes, k=factors, scIdx=groupings, isLog=FALSE)


write.table(ruv$normalizedCounts, "matrix.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=F)
write.table(ruv$W, "ruv_weights.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=F)
write.table(groupings, "makeGroups.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=F)

#plotting
colors <- rep("gray", length(factorX))
colors[which(factorX=="replicate_control")] = "red"
pdf("Rplots.pdf")
plotRLE(matx, outline=FALSE, ylim=c(-4,4), las=2, main="Pre RUV: Relative Log Expression Plot",col= colors)
plotPCA(matx, cex=1.2, col=colors, main="Pre RUV: PCA Plot")
plotRLE(ruv$normalizedCounts, outline=FALSE, ylim=c(-4,4), las=2, main="Post RUV: Relative Log Expression Plot", col=colors)
plotPCA(matx, cex=1.2, col=colors, main="Post RUV: PCA Plot")
dev.off()

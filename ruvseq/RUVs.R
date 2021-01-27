#!/usr/bin/env R

usage = "Rcript /Users/patterja/galaxy/tools/compbio-galaxy-wrappers/ruvseq/RUVs.R /Users/patterja/galaxy/database/files/000/dataset_215.dat FALSE 3 3 1 all"

version="0.2.0"

ruvSeqPackageExists <- require ("RUVSeq")
if ( !ruvSeqPackageExists ) {install.packages ("RUVSeq")}

library(RUVSeq)

#Get the inputs
args <- commandArgs(TRUE)
dataFile <- args[1]
normalize = args[2]
minNumReads <- as.double(args[3])
minNumSamples <- as.double(args[4])
idxSamples <-args[5]
factors <- as.double(args[6])
repID <- args[7]
cIdx <- args[8]

#read datafile
data <- read.table(dataFile, header = TRUE, sep="\t", row.names=1, check.names=FALSE)

#BCCL is included then use counts and normalize
if (normalize==TRUE){
    cpm = apply(data,2, function(x) (x/sum(x))*1000000)
    ndata = cpm
} else {
    ndata = data
}

#filter for smmart samples that are larger then thresholds
# TODO: filter that only includes UHR only/samples only (not both)
# i don't know of a better way to do this. 


smmartsamps = colnames(ndata)[eval(parse(text=idxSamples))]
smmartcpm = ndata[,smmartsamps]
#FILTER
filter <-  apply(smmartcpm, 1, function(x) length(x[x>minNumReads])>=minNumSamples)
filtered <- ndata[filter,]

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


#prepare for ruvs
matx <- as.matrix(filtered)
matfiltered <- round(matx, 0)

#identify replicates for ruvs input
samp_reps <- colnames(filtered)
reps = unlist(lapply(unlist(strsplit(repID,",")), trimws))
for (rep in reps){
   samp_reps[grep(rep, samp_reps, value=FALSE)] = rep
}

factorX <- as.factor(c(samp_reps))
groupings <- makeGroups(factorX)

ruv <- RUVs(matfiltered, cIdx=genes, k=factors, scIdx=groupings, isLog=FALSE)


#output
write.table(ruv$normalizedCounts, "normalized.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=F)
write.table(ruv$W, "ruv_weights", sep="\t", row.names=TRUE, col.names=TRUE, quote=F)
write.table(groupings, "groupings", sep="\t", row.names=FALSE, col.names=FALSE, quote=F)
#plotting
colors <- rep("gray", length(factorX))
colors[which(factorX=="replicate_control")] = "red"
pdf("Rplots.pdf")
plotRLE(matx, outline=FALSE, ylim=c(-4,4), las=2, main="Pre RUV: Relative Log Expression Plot",col= colors)
plotPCA(matx, cex=1.2, col=colors, main="Pre RUV: PCA Plot")
plotRLE(ruv$normalizedCounts, outline=FALSE, ylim=c(-4,4), las=2, main="Post RUV: Relative Log Expression Plot", col=colors)
plotPCA(matx, cex=1.2, col=colors, main="Post RUV: PCA Plot")
dev.off()

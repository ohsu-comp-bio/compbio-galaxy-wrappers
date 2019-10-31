

#### Package Install ####
ruvSeqPackageExists <- require ("RUVSeq")

if ( !ruvSeqPackageExists ) {
  install.packages ("RUVSeq")
  library ("RUVSeq")
  
}

# rcolorPackageExists <- require ("RColorBrewer")

# if ( !rcolorPackageExists ) {
  # install.packages ("RColorBrewer")
  # library ("RColorBrewer")
  
# }

#########################

#Get the inputs
args <- commandArgs(TRUE)
dataFile <- args[1]
factors <- as.double(args[2])
minNumReads <- as.double(args[3])
minNumSamples <- as.double(args[4])
repID <- args[5]
cIdx <- args[6]


phenoData <- read.table(dataFile, header = TRUE, sep="\t")

filter <-  apply(phenoData, 1, function(x) length(x[x>minNumReads])>=minNumSamples)
filtered <- phenoData[filter,]

rownames(filtered) <- filtered[,1]

filtered <- filtered[-c(1)]

genes <- rownames(filtered)[grep(cIdx, rownames(filtered))]

if (cIdx=="all"){
	genes <- rownames(filtered)
}

if (grepl(",", cIdx, fixed=TRUE)) {
	genes <- strsplit(cIdx, ",")
	genes <- unlist(genes, use.names=FALSE)
}

x <- colnames(filtered)

for (index in 1:length(x)) {
	name=x[index]
	if (regexpr(repID, name)!=-1){
		x[index]=repID
	}
}


matrix <- as.matrix(filtered)

counts <- round(matrix, 0)
set <- newSeqExpressionSet(counts, phenoData = data.frame(x, row.names=colnames(filtered)))

factorX <- as.factor(c(x))
differences <- makeGroups(factorX)

set1 <- RUVg(counts, genes, k=factors)
set2 <- RUVg(set, genes, k=factors)

capture.output(pData(set2), file="pData.txt")

write.table(set1$normalizedCounts, "matrix.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=F)
write.table(set1$W, "ruv_weights.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=F)
write.table(differences, "makeGroups.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=F)


#colors <- brewer.pal(3, "Set2")
#plotRLE(set2, outline=FALSE, ylim=c(-4,4), col=colors[x])
#plotPCA(set2, col=colors[x], cex=1.2)

plotRLE(set2, outline=FALSE, ylim=c(-4,4), las=2, main="Relative Log Expression Plot")
plotPCA(set2, cex=1.2, main="PCA Plot")

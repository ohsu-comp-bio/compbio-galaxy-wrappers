
args <- commandArgs(TRUE)
target_id <- args[1]
centroids <- args[2]
input <- args[3]
hugo <- as.logical(args[4])

##------------------GENE NAME CONVERTER FILE. It converts ENSG names to HUGO gene names 
targetid = read.csv(target_id, sep="\t", stringsAsFactors = F)
targetid <- unique(targetid[, c("X.1", "X.5")])
rownames(targetid) <- targetid$X.1

##------------------pam50 centroid file
pam50_centroids = read.csv(centroids, check.names = FALSE, sep="\t", row.names = 1)


##------------------Matrix file: smmart samples matrix (columns samples, rows genes)

smmart_tpm = read.csv(input, sep="\t", row.names = 1, header=TRUE)

##------------------CONVERT ENSG NAMES TO HUGO GENE NAMES
if ( !hugo ){
	smmart_tpm = smmart_tpm[intersect(rownames(smmart_tpm), targetid$X.1),]
	smmart_tpm <- as.matrix(smmart_tpm)
}
rownames(smmart_tpm) = targetid[intersect(rownames(smmart_tpm), targetid$X.1), "X.5"]


#there are a few ways to sum the values based on rownames. There are probably faster functions for this
smmart_tpm = aggregate(smmart_tpm,by=list(rownames(smmart_tpm)), FUN=sum)


#turn the first column into the row names
rownames(smmart_tpm) <- smmart_tpm$Group.1
smmart_tpm <- smmart_tpm[,-1]

##------------------LOG2+1
# log2(x+1) the ruv fixed matrix
l.smmartexpmat=log2(smmart_tpm+1)

##------------------Z-SCALING: Calculate z-score for each gene across all samples: z-score=(value-gene mean)/gene sd
sctl_smmartexpmat = t(scale(data.frame(t(l.smmartexpmat), check.names =F), scale=T))

#samples with sd=0 will calc a z-score of NA (dividing by zero) remove these rows
sctl_smmartexpmat = sctl_smmartexpmat[complete.cases(sctl_smmartexpmat),]


#turn to matrix
mat_tpm = as.matrix((sctl_smmartexpmat))



## CORRELATION WITH CENTROIDS
scor=round(cor(mat_tpm[intersect(rownames(mat_tpm), rownames(pam50_centroids)),], pam50_centroids[intersect(rownames(mat_tpm), rownames(pam50_centroids)),], method="pearson"),2)

#write to output
write.table(smmart_tpm, file="smmart_tmp.tsv", row.names=T, col.names=NA, sep="\t", quote=F)

write.table(scor, file="scor.tsv", row.names=T, sep="\t", col.names=NA, quote=F)


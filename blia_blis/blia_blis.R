# VERSION = 0.1.0

library(ggplot2)
library(reshape2)
library(stats)
library(psych)
library(RUVSeq)

plot_me <- function (name, to_plot) {
  pdf(name)
  plot(to_plot)
  dev.off()
}

# Set arguments.
# blia_blis.R <meta> <df_hugo> <centroids> <cIdx> <batch>
args <- commandArgs(trailingOnly=TRUE)

# meta <- read.csv(args[1], header=TRUE, sep="\t")
dataFile <- args[1]
df_hugo_raw <- read.csv(dataFile, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
# centroids <- read.csv(args[3], header=TRUE, row.names=1, sep="\t")
cIdx <- read.csv(args[2], header=FALSE)
# batch <- read.csv(args[5], header=TRUE, sep="\t")

# Correlation against passed centroids, create plot with CIs.
# df_hugo <- df_hugo_raw[,as.character(meta$Accession)]
# colnames(df_hugo) <- paste(meta$Patient, meta$Biopsy, sep="_")
# tnbc_res <- cor(df_hugo[rownames(centroids),], centroids, method="spearman")
# res <- corr.test(as.matrix(df_hugo[rownames(centroids),]), as.matrix(centroids), method="spearman")
# m1 <- do.call(rbind, strsplit(rownames(res$ci), '-'))
# res$ci$Sample <- paste(m1[,1], m1[,2], sep="-")
# res$ci$Subtype <- m1[,3]
# ggplot(res$ci, aes(x=Sample, y=r, fill=Subtype,label = ifelse(p < 0.01, "*", ""))) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_errorbar( aes(ymin=lower, ymax=upper), inherit.aes = TRUE, position = "dodge") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90)) +
#   guides(fill=guide_legend(ncol=1)) +
#   geom_text(vjust = -7, position = position_dodge(width = 1))

factors <- 2
minNumReads <- 3
minNumSamples <- 3
repID <- 'UHR'
phenoData <- read.table(dataFile, header = TRUE, sep="\t", row.names=1, check.names=FALSE)
filter <-  apply(phenoData, 1, function(x) length(x[x>minNumReads])>=minNumSamples)
filtered <- phenoData[filter,]
genes <- strsplit(cIdx[,1], ",")
genes <- unlist(genes, use.names=FALSE)

matx <- as.matrix(filtered)
matfiltered <- round(matx, 0)

samp_reps <- colnames(filtered)
reps <- unlist(lapply(unlist(strsplit(repID,",")), trimws))
for (rep in reps){
  samp_reps[grep(rep, samp_reps, value=FALSE)] <- rep
}

factorX <- as.factor(c(samp_reps))
groupings <- makeGroups(factorX)

ruv <- RUVg(matfiltered, cIdx=genes, k=factors, isLog=FALSE)

colors <- rep("gray", length(factorX))
colors[which(factorX=="replicate_control")] <- "red"

plot_batch_tpm <- function(au, batch){
  au <- log(au+1, base=2)
  au <- au[,1:ncol(au)]
  au <- reshape2::melt(as.matrix(au))
  colnames(au) <- c("Gene", "Sample", "TPM")
  au$Sample <- as.character(au$Sample)
  au <- merge(au, batch, by="Sample")
  ggplot(au, ggplot2::aes(x=Sample, y=TPM, color=Batch)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.text = element_text(size = 4), legend.position = "bottom", legend.title = element_text(size = 6)) +
    ylab("Log2 TPM")
}

au_un <- df_hugo_raw
au <- ruv$normalizedCounts
# filtered by read count
filter <-  apply(au, 1, function(x) length(x[x>minNumReads])>=minNumSamples)
au2 <- au[filter,]
filter <-  apply(au_un, 1, function(x) length(x[x>minNumReads])>=minNumSamples)
au2_un <- au_un[filter,]

# # Create the log2 TPM plots
# plot_batch_tpm(au_un, batch)
# plot_batch_tpm(au2_un, batch)
# plot_batch_tpm(au, batch)
# plot_batch_tpm(au2, batch)

# This stuff is going to Rplots since I don't know how to write them to file
par(mar=c(7,4,4,2)+.1)
plotRLE(matx, outline=FALSE, cex.axis=0.9, ylim=c(-4,4), xaxt="n", main="Pre RUV: Relative Log Expression Plot")
labs <- colnames(matx)
axiscolors <- c(rep("black", 33), "blue")
Map(axis, side=1, at=1:34, cex.axis=0.6, font=2, col.axis=axiscolors, labels=labs, lwd=0, las=2)
axis(1, at=1:34, labels=FALSE)

plotRLE(ruv$normalizedCounts, outline=FALSE, cex.axis=0.9, ylim=c(-4,4), xaxt="n", main="Pre RUV: Relative Log Expression Plot")
labs <- colnames(ruv$normalizedCounts)
axiscolors <- c(rep("black", 33), "blue")
Map(axis, side=1, at=1:34, cex.axis=0.6, font=2, col.axis=axiscolors, labels=labs, lwd=0, las=2)
axis(1, at=1:34, labels=FALSE)

par(mar=c(5,5,4,2)+.1)
plotPCA(matx, cex=0.6, main="Pre RUV: PCA Plot")
plotPCA(ruv$normalizedCounts, cex=0.6, main="Post RUV: PCA Plot")

# Text output.
write.table(ruv$normalizedCounts, "normalized.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=F)
write.table(matx, "unnormalized.tsv", sep="\t", row.names=TRUE, col.names=NA, quote=F)
write.table(ruv$W, "ruv_weights", sep="\t", row.names=TRUE, col.names=TRUE, quote=F)
write.table(groupings, "groupings", sep="\t", row.names=FALSE, col.names=FALSE, quote=F)

library(corrplot)
library(pheatmap)
options(warn=-1)

dg_corr_plot <- function(X){
  row_order <- sort(rowSums(X^2, na.rm = TRUE), decreasing = TRUE, index.return = TRUE)$ix
  col_order <- sort(colSums(X^2, na.rm = TRUE), decreasing = TRUE, index.return = TRUE)$ix
  X <- X[row_order, col_order, drop = FALSE]
  X <- X[which(rowMeans(is.na(X)) != 1),]
  X[is.na(X)] <- 0
  corrplot(X, tl.col = "black", tl.cex = .50, method = "circle", cl.pos = "n", col = corr_colors)
}

dg_corr_table <- function(X){
  X <- X[which(rowMeans(is.na(X)) != 1),]
  X <- cbind.data.frame(as.vector(matrix(rownames(X), nrow = nrow(X), ncol = ncol(X), byrow = FALSE)), as.vector(matrix(colnames(X), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)), as.vector(X))
  colnames(X) <- c("Drug", "Gene", "Correlation")
  X[is.na(X)] <- 0
  X <- X[sort(abs(X[,"Correlation"]), decreasing = TRUE, index.return = TRUE, na.last = NA)$ix,]
  return (X)
}

dd_corr_plot <- function(X){
  X[is.na(X)] <- 0
  corrplot(X, tl.col = "black", tl.cex = .35, method = "circle", cl.pos = "n", order = "hclust", col = corr_colors)
}

dd_corr_table <- function(X){
  X[lower.tri(X)] <- NA
  X <- cbind.data.frame(as.vector(matrix(rownames(X), nrow = nrow(X), ncol = ncol(X), byrow = FALSE)), as.vector(matrix(colnames(X), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)), as.vector(X))
  colnames(X) <- c("Drug #1", "Drug #2", "Correlation")
  X <- X[which(!is.na(X[, 3])), ]
  X <- X[which(X[,1] != X[,2]),]
  X <- X[sort(abs(X[,"Correlation"]), decreasing = TRUE, index.return = TRUE)$ix,]
  return (X)
}

args <- commandArgs(TRUE)
# set up data
drug_screen <- read.csv(args[1], stringsAsFactors = FALSE, sep=",", check.names = FALSE)
# drug_screen2 <- read.csv(args[2], stringsAsFactors = FALSE, sep=",", check.names = FALSE)
# drug_screen3 <- read.csv(args[3], stringsAsFactors = FALSE, sep=",", check.names = FALSE)
# drug_screen_tmp <- rbind(drug_screen1, drug_screen2)
# drug_screen <- rbind(drug_screen_tmp, drug_screen3)

RNASeq <- read.csv(args[2], stringsAsFactors = FALSE, sep=",", row.names = "gene", check.names = FALSE)
# RNASeq2 <- read.csv(args[5], stringsAsFactors = FALSE, sep=",", row.names = "gene", check.names = FALSE)
# RNASeq3 <- read.csv(args[6], stringsAsFactors = FALSE, sep=",", row.names = "gene", check.names = FALSE)
# RNASeq_tmp <- cbind(RNASeq1, RNASeq2)
# RNASeq <- cbind(RNASeq_tmp, RNASeq3)

AUC_All <- tapply(drug_screen[,"Inhibitor Interpreted Result: Area under the curve"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)
IC50_All <- tapply(drug_screen[,"Inhibitor Interpreted Result: IC50"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)

BeatAML_DD_All <- list()
BeatAML_DD_All$AUC <- AUC_All
BeatAML_DD_All$IC50 <- IC50_All
BeatAML_DD_All$AUC_AUC_Pearson <- cor(AUC_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_All$AUC_AUC_Spearman <- cor(AUC_All, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DD_All$IC50_IC50_Pearson <- cor(IC50_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_All$IC50_IC50_Spearman <- cor(IC50_All, method = "spearman", use = "pairwise.complete.obs")

RNASeq <- log2(t(RNASeq) + 1)
common_patients_All <- intersect(rownames(IC50_All), rownames(RNASeq))

gene_names <- c("CTLA-4", "Galectin-9", "PD-L1", "PD-L2", "PD1", "TIM-3", "VISTA", "LAG3", "TNFRSF4", "CD3G", "ARG1", "ARG2")
gene_ids <- c("CTLA4", "LGALS9", "CD274", "PDCD1LG2", "PDCD1", "HAVCR2", "C10orf54", "LAG3", "TNFRSF4", "CD3G", "ARG1", "ARG2")
RNASeq_All <- RNASeq[common_patients_All, gene_ids]
colnames(RNASeq_All) <- gene_names

AUC_All <- AUC_All[common_patients_All,]
IC50_All <- IC50_All[common_patients_All,]

BeatAML_DG_All <- list()
BeatAML_DG_All$AUC <- AUC_All
BeatAML_DG_All$IC50 <- IC50_All
BeatAML_DG_All$RNASeq <- RNASeq_All
BeatAML_DG_All$AUC_RNASeq_Pearson <- cor(AUC_All, RNASeq_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_All$AUC_RNASeq_Spearman <- cor(AUC_All, RNASeq_All, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DG_All$IC50_RNASeq_Pearson <- cor(IC50_All, RNASeq_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_All$IC50_RNASeq_Spearman <- cor(IC50_All, RNASeq_All, method = "spearman", use = "pairwise.complete.obs")

input_genes <- c("CTLA-4", "Galectin-9", "PD-L1", "PD-L2", "PD1", "TIM-3", "VISTA", "LAG3", "TNFRSF4", "CD3G", "ARG1", "ARG2")

colors <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
corr_colors <- rev(colors(201))

if(args[3] == 'pearson' && args[4] == 'AUC'){
  
  pdf(args[5], width=7, height=10.5)
  dg_corr_plot(BeatAML_DG_All$AUC_RNASeq_Pearson[, input_genes, drop = FALSE])
  dev.off()
  
  pdf(args[6], width=12, height=8, paper="A4r")
  dd_corr_plot(BeatAML_DD_All$AUC_AUC_Pearson)
  dev.off()
  
}else if(args[3] == 'spearman' && args[4] == 'AUC'){
  
  pdf(args[5], width=7, height=10.5)
  dg_corr_plot(BeatAML_DG_All$AUC_RNASeq_Spearman[, input_genes, drop = FALSE])
  dev.off()
  
  pdf(args[6], width=12, height=8, paper="A4r")
  dd_corr_plot(BeatAML_DD_All$AUC_AUC_Spearman)
  dev.off()
  
}else if(args[3] == 'pearson' && args[4] == 'IC50'){
  
  pdf(args[5], width=7, height=10.5)
  dg_corr_plot(BeatAML_DG_All$IC50_RNASeq_Pearson[, input_genes, drop = FALSE])
  dev.off()
  
  pdf(args[6], width=12, height=8, paper="A4r")
  dd_corr_plot(BeatAML_DD_All$IC50_IC50_Pearson)
  dev.off()
  
}else if(args[3] == 'spearman' && args[4] == 'IC50'){
  
  pdf(args[5], width=7, height=10.5)
  dg_corr_plot(BeatAML_DG_All$IC50_RNASeq_Spearman[, input_genes, drop = FALSE])
  dev.off()
  
  pdf(args[6], width=12, height=8, paper="A4r")
  dd_corr_plot(BeatAML_DD_All$IC50_IC50_Spearman)
  dev.off()
}


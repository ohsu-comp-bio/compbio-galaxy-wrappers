#!/usr/bin/env Rscript

# Current Version: 0.1.1
# Version history
# 0.1.0 - initial commit
# 0.1.1 - set output file names as args

### Load required packages

library(argparser)
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(edgeR)
library(reshape)
library(rrcov)
library(ggfortify)
library(ggthemes)
library(gridExtra)


### Collect command line arguments

# create parser
p <- arg_parser("UHR Quality Control")

# add command line arguments
p <- add_argument(p, "background_cohort", help="text file defining background cohort samples (.txt)")
p <- add_argument(p, "validation_targets", help="file defining list of validation targets for cancer related genes (.tsv)")
p <- add_argument(p, "new_sample_run", help="sequencing run ID of new sample (string)")
p <- add_argument(p, "new_sample_counts", help="kallisto output for new sample (.tsv)")
p <- add_argument(p, "gene_set", help="list of cancer related genes of interest (.txt)")
p <- add_argument(p, "target_map", help="file that maps genes to target_ids (.tsv)")
p <- add_argument(p, "counts_matrix", help="Kallisto counts matrix with all UHR samples (.tsv)")
p <- add_argument(p, "qc_report", help="file name for plots output (.pdf)")
p <- add_argument(p, "stats_report", help="file name for statistics output (.txt)")

# parse the command line arguments
argv <- parse_args(p)

# read data supplied by arguments
background_samples <- scan(argv$background_cohort, what="", sep="\n")

valid_targets <- read_tsv(argv$validation_targets, show_col_types=F)["target_id"]
valid_targets <- as.list(valid_targets$target_id)

new_run <- argv$new_sample_run

new_sample <- read_tsv(argv$new_sample_counts, show_col_types=F)["est_counts"]
colnames(new_sample) <- new_run

gene_set <- as.list(scan(argv$gene_set, what="", sep="\n"))

target_map <- read_tsv(argv$target_map, show_col_types=F)
colnames(target_map) <- c("target_id","gene")

raw_counts <- read_tsv(argv$counts_matrix, show_col_types=F)

qc_report <- argv$qc_report
stats_report <- argv$stats_report

### Build matrix of raw Kallisto counts data containing background cohort samples and new sample

# subset full Kallisto counts matrix to background cohort samples
coh_raw <- raw_counts[,background_samples]

# get genomic coordinates for Kallisto data
target_id <- colsplit(raw_counts$target_id, split = "\\|", names = c("target_id"))[1]

# add target coordinate and new sample to finalize raw Kallisto counts matrix
raw_counts <- cbind(target_id, coh_raw, new_sample)


### Subset raw data to cancer genes of interest

# make map of cancer specific genes to targets
cancer_targs <- data.frame(matrix(ncol = 1, nrow = 0))
for(name in gene_set){

  # find targets associated with specific gene name
  tmp_targets <- subset(target_map, gene == name)
  cancer_targs <- rbind(cancer_targs, tmp_targets)
}

# subset counts data to cancer genes of interest
sub_df <- raw_counts[raw_counts[["target_id"]] %in% cancer_targs[["target_id"]],]

# add gene names to dataframe
sub_raw <- data.frame(matrix(ncol = 1, nrow = 0))
for (name in gene_set){

  # find targets associated with specific gene name
  tmp_targets <- subset(target_map, gene == name)
  tmp_df <- sub_df[sub_df[["target_id"]] %in% tmp_targets[["target_id"]],]

  # add gene name
  if(length(tmp_df != 0)){gene = name; tmp_df <- cbind(gene, tmp_df)}

  # append to dataframe
  sub_raw <- rbind(sub_raw,tmp_df)
}


### TMM normalization of raw counts data

# build DGE counts object
dgeCounts <- DGEList(subset(sub_raw, select = -c(gene, target_id)), genes = sub_raw["target_id"])

# calculate normalization factors using the TMM method
tmm_factors <- calcNormFactors(dgeCounts, method="TMM")

# create dataframe of normalized counts using TMM normalization
sub_tmm <- data.frame(cpm(tmm_factors, log=T), check.names = FALSE)
sub_tmm <- cbind(sub_raw["gene"], sub_raw["target_id"], sub_tmm)


### Compare new sample to background cohort

# separate dataframe into background cohort and new sample
background_tmm <- sub_tmm[-c(ncol(sub_tmm))]
new_sample_tmm <- sub_tmm[c("gene","target_id", new_run)]


# subset each dataframe to gene validation targets
new_sample_valid <- new_sample_tmm[new_sample_tmm[["target_id"]] %in% valid_targets,]
new_sample_valid_data <- subset(new_sample_valid, select = -c(target_id))
background_valid <- background_tmm[background_tmm[["target_id"]] %in% valid_targets,]
background_valid_data <- subset(background_valid, select = -c(target_id))

# initialize variables to be used in for loop
ci_test <- logical(0)
z_scores <- numeric(0)
plot_lst <- list()

# build plots for new point against background distribution for all gene validation targets
validation_genes <- as.list(background_valid$gene)

melted <- (matrix(ncol = 0, nrow = 1))
for(name in validation_genes){

  # subset data to specified gene
  tmp <- subset(background_valid_data, gene == name)
  tmp_new <- subset(new_sample_valid_data, gene == name)
  melted <- melt(tmp, id = "gene")

  # calculate 95% CI
  median <- rowMedians(as.matrix(subset(tmp, select = -c(gene))))
  sd <- rowSds(as.matrix(subset(tmp, select = -c(gene))))
  lower <- median - 1.96 * sd
  upper <- median + 1.96 * sd

  # test if new data point is within the 95% CI
  inside = F
  if(tmp_new[new_run] >= lower & tmp_new[new_run] <= upper){inside = T}
  ci_test <- c(ci_test, inside)

  # calculate z-score of new data point against background distribution
  score <- as.numeric((tmp_new[new_run]-median)/sd)
  z_scores <- c(z_scores, score)

  # plot distribution of each validation target and how new sample aligns
  p <- ggplot(melted, aes(x = value)) +
    ggtitle(name) +
    geom_histogram(aes(y=after_stat(density)), bins = 12, color = "black", fill = "grey84") +
    xlab("Counts per million (log2)") +
    geom_density(color = "red") +
    geom_vline(aes(xintercept = as.numeric(tmp_new[ncol(tmp_new)])) , tmp_new, color = "blue", linewidth = 2) +
    annotate("rect", xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf, fill = "darkgrey", alpha = 0.3) +
    theme_clean()

  # add plot to list of ggplots
  plot_lst[[name]] <- ggplotGrob(p)

}

# build file to store text outputs
sink(stats_report)

# calculate principle components
trans <- data.frame(t(subset(sub_tmm, select = -c(gene, target_id))))
sample <- data.frame(as.factor(rownames(trans)))
colnames(sample) <- "sample"
sample["group"] <- "background cohort"
sample$group[length(sample$group)] <- "new sample"
pca <- prcomp(trans, scale = T)

# calculate Spearman correlation coefficient between new sample and background cohort median values
background <- sub_tmm[-c(1,2,ncol(sub_tmm))]
background["median"] <- rowMedians(as.matrix(background))
cor.test(background[["median"]], sub_tmm[[new_run]], method = "spearman")

# calculate pairwise Spearman correlation coefficients between all UHRs
cor_matrix <- cor(as.matrix(subset(sub_tmm, select = -c(gene, target_id))), method = "spearman")

# calculate Pearson correlation coefficient between new sample and background cohort median values
dat <- data.frame(background[["median"]], sub_tmm[[new_run]])
colnames(dat) <- c("median", "new")
cor.test(background[["median"]], sub_tmm[[new_run]])

# proportion of validation targets in the 95% CI (median +/- 1.96 * SD)
prop <- round(sum(ci_test)/length(ci_test),3)
cat("\nProportion of new sample validation targets within the 95% CI: ", prop)

# calculate mean and sd of z-scores
z_scores <- data.frame(z_scores)
z_mean <- round(mean(z_scores$z_scores),3)
z_sd <- round(abs(sd(z_scores$z_scores)),3)
cat("\nZ-score Mean: ", z_mean)
cat("\nZ-score SD", z_sd, "\n")
sink()

### Build PDF file with graphics and outputs

pdf(file = qc_report, width = 10, height = 7)

# add PCA plot to pdf
autoplot(pca, data = sample, colour = "group", main = "PCA")

# add pairwise Spearman correlation heatmap to pdf
melted_cor_mat <- melt(cor_matrix)
ggplot(data = melted_cor_mat, aes(x=X1, y=X2,fill=value)) +
  geom_tile() +
  coord_fixed() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("Pairwise Spearman Correlation Heatmap")

# add correlation scatterplot to pdf
ggplot(dat, aes(median, new)) +
  ggtitle("Correlation Scatterplot") +
  xlab("Background Cohort Median") +
  ylab("New Sample") +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

# plot distribution of z-scores for new sample against background cohort for Rodrigo validation targets
header <- "Distribution of Validation Target Z-Scores for New Sample"
ggplot(z_scores, aes(z_scores)) +
  ggtitle(header) +
  geom_histogram(aes(y=after_stat(density)),bins=20, color = "black", fill = "grey84") + xlab("Z Score") +
  geom_density(color = "red")

# make multiple boxplot with distribution of background cohort and new sample comparison
melted = melt(background_valid_data)
ggplot(melted, aes(x = gene, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = new_sample_valid, aes(x = gene, y = new_sample_valid[[new_run]], color = "red"), show.legend = F) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_clean() +
  ggtitle("New Sample Against Background Cohort") +
  xlab("Gene") +
  ylab("Counts per Million (log2)")

# add distribution plots (histogram) of validation targets to pdf in grid format (3x3)
marrangeGrob(plot_lst, ncol = 3, nrow = 3)

dev.off()


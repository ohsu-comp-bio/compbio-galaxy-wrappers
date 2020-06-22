#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser()
parser$add_argument("samples", help="File with sample names")
parser$add_argument("data_dir", help="Input directory with VCF files for samples")
parser$add_argument("signatures", help="File of cosmic signature probabilities")
parser$add_argument("--plotbysample", action="store_true", default=FALSE, help="Whether to plot individual")
parser$add_argument("--sig3", default=NULL, help="Threshold to set for Signature 3 Contribution")
args <- parser$parse_args()


# Load ref
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Load samples
sample_names <- scan(file=args$samples, what=character())
print(sample_names)
vcf_files <- unlist(sapply(sample_names, function(x){list.files(args$data_dir, pattern = x, full.names = TRUE)}))
print(vcf_files)
vcfs <- read_vcfs_as_granges(vcf_files, names(vcf_files), ref_genome)

# Plot Mutation Spectrum
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
pdf("cosmic_signatures_type_occurences.pdf", width=7,height=7)
plot_spectrum(type_occurrences, CT = TRUE)
dev.off()


# Mutation Profile
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

if (args$plotbysample == TRUE) {
	pdf("cosmic_signatures_type_occurences_by_sample.pdf", width=7,height=7)
	print(plot_spectrum(type_occurrences, by = sample_names, CT = TRUE, legend = TRUE))
	dev.off()

	pdf("cosmic_signatures_profiles.pdf", width=7,height=7)
	print(plot_96_profile(mut_mat))
	dev.off()
}


# Cosmic Signatures
cancer_signatures = read.table(args$signatures, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])
pdf("cosmic_signatures_profiles_by_signature.pdf", width=7,height=7)
plot_96_profile(cancer_signatures[,1:3], condensed = TRUE, ymax = 0.3)
dev.off()

# Fit to cosmic signatures
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(fit_res$contribution) > 10)

# Write results
write.table(fit_res$contribution,file="cosmic_sig_fit_contributions.tsv", sep="\t", quote=FALSE, col.names=NA)

pdf("cosmic_signatures_fit_barplot.pdf", width=7,height=7)
plot_contribution(fit_res$contribution[select,,drop=FALSE],cancer_signatures[,select,drop=FALSE], coord_flip = FALSE,mode = "absolute")
dev.off()

if (length(vcf_files) > 4) {
    pdf("cosmic_signatures_fit_heatmap.pdf", width=8, height=10)
    plot_contribution_heatmap(fit_res$contribution, cluster_samples = TRUE,method = "complete")
    dev.off()
}

# Filter and write signature 3 samples and contributions
if (!is.null(args$sig3)) {
	maxs <- apply(a, 2, max)
	mins <- apply(a, 2, min)
	#relative <- scale(fit_res$contribution, center = mins, scale = maxs - mins)
	sig3 <- fit_res$contribution[3,fit_res$contribution$Signature.3 > args$sig3, drop=FALSE]
	write.table(sig3,file="cosmic_sig3_contributions.tsv", sep="\t", quote=FALSE, col.names=NA)
}

# Make sample barplots
if (args$plotbysample == TRUE) {
	for (sample in 1:length(sample_names)){
		samp = sample_names[sample]
		contri <- fit_res$contribution[,samp, drop=FALSE]
		data <- melt(as.matrix(contri))
		colnames(data) <- c("Signature", "Sample", "Contribution")
		g <- ggplot(data, aes(x=Signature, y=Contribution)) +
		    geom_bar(stat="identity", position='dodge', fill="#2b8cbe") +
		    theme_minimal() +
		    theme(axis.text.x = element_text(angle = 90, hjust = 1))

		ggsave(paste("cosmic_signatures_",samp,"_barplot.png",sep=""),g,width=7,height=7)
	}
}




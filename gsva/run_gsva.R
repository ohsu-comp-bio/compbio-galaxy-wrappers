suppressWarnings(suppressMessages(library(GSVA)))
suppressWarnings(suppressMessages(library(GSEABase)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(dplyr)))

# Wrapper for GSVA

# TODO (hocumj): Implement proper error and warning handling.
# TODO (hocumj): Allow the user to set all optional arguments.

args <- commandArgs(TRUE)

input_file <- args[1]
gene_set_filename <- args[2]
min_size_set <- as.integer(args[3])
max_size_set <- as.integer(args[4])
expression_cutoff_quantifications <- as.numeric(args[5])
expression_cutoff_samples <- as.integer(args[6])
log_transform_expression <- as.integer(args[7])
kernel <- as.character(args[8])
max_diff <- as.integer(args[9])
num_of_cores <- as.integer(args[10])
verbose <- as.integer(args[11])
output_file <- args[12]

if (log_transform_expression == 1)
{
    log_transform_expression <- TRUE
} else 
{
    log_transform_expression <- FALSE
}
if (max_diff == 1)
{
    max_diff <- TRUE
} else 
{
    max_diff <- FALSE
}
if (verbose == 1)
{
    verbose <- TRUE
} else 
{
    verbose <- FALSE
}

gene_expression_values = read.csv(input_file, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
# Prepare gene sets.
gene_sets <- GeneSetCollection(c(getGmt(gene_set_filename)))

# Subset genes based on expression.
# Threshold for minimum expression and minimum number of samples expression threshold must be met.
gene_expression_values <- gene_expression_values[apply(gene_expression_values[, -1], 1, function(X) length(X[X >= expression_cutoff_samples]) > expression_cutoff_quantifications), ]

if (log_transform_expression == TRUE)
{
    gene_expression_values <- log(gene_expression_values, base = 2)
}

# Calculate Enrichment Scores with GSVA.
# TODO (hocumj): Allow passing all optional arguments such as tau and method.
gsva_results <- suppressWarnings(gsva(as.matrix(gene_expression_values), gene_sets, min.sz = min_size_set, max.sz = max_size_set, kcdf = kernel, mx.diff = max_diff, verbose = verbose, parallel.sz = num_of_cores))

# Prepare GSVA output.
gsva_results <- melt(gsva_results)
colnames(gsva_results) <- c("Pathway", "Sample", "EnrichmentScore")

write.table(gsva_results, file = output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
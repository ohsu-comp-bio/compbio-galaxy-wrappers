# VERSION = 2.0.0

jsonExists <- require("jsonlite")

if (!jsonExists){
	install.packages("jsonlite")
	library("jsonlite")
}

args <- commandArgs(TRUE)

file <- args[1]
sample_name <- args[2]

sample_data <- fromJSON(file)

print(sample_data)

validation_data <- data.frame("TPM" = c(0.01, 0.1, 1, 10, 100, 1000), "Num_genes" = c(17354, 15828, 12618, 6528, 741, 58), "sdev" = c(608, 647, 673, 1093, 241, 10))

chart_data <- validation_data

chart_data$sample_tpm <- 0

print(chart_data)

sample_data <- sample_data$sampleRunMetrics

print(sample_data)

print(sample_data[sample_data$metric == "rna_tpm_zero", c("value")])

chart_data[chart_data$TPM == 0.01, c("sample_tpm")] <- sample_data[sample_data$metric == "rna_tpm_hundredth", c("value")]
chart_data[chart_data$TPM == 0.1, c("sample_tpm")] <- sample_data[sample_data$metric == "rna_tpm_tenth", c("value")]
chart_data[chart_data$TPM == 1, c("sample_tpm")] <- sample_data[sample_data$metric == "rna_tpm_one", c("value")]
chart_data[chart_data$TPM == 10, c("sample_tpm")] <- sample_data[sample_data$metric == "rna_tpm_ten", c("value")]
chart_data[chart_data$TPM == 100, c("sample_tpm")] <- sample_data[sample_data$metric == "rna_tpm_onehundred", c("value")]
chart_data[chart_data$TPM == 1000, c("sample_tpm")] <- sample_data[sample_data$metric == "rna_tpm_onethousand", c("value")]

#chart_data$TPM <- format(chart_data$TPM, scientific = FALSE)

print(chart_data)

ymax = as.numeric(max(chart_data$Num_genes[1], chart_data$sample_tpm[1]))

plot(chart_data$TPM, chart_data$Num_genes, type="n", xlab="Expression threshold (log10 TPM)", ylab="Number of genes", log="x", xaxt="n", tck=1, ylim=c(0,ymax))

axis(1, at=c(0.01, 0.1, 1, 10, 100, 1000), labels=c("0.01", "0.1", "1", "10", "100", "1000"), tck=1)



lines(chart_data$TPM, chart_data$Num_genes, type="b", pch=15, lty=3, col="blue", lwd=2)
arrows(chart_data$TPM, (chart_data$Num_genes - chart_data$sdev), chart_data$TPM, (chart_data$Num_genes + chart_data$sdev), length=0.05, angle=90, code=3, lwd=1.5)
lines(chart_data$TPM, chart_data$sample_tpm, type="b", pch=19, lty=1, col="orange", lwd=2)


legend("topright", legend=c("Validation Sample Average", sample_name), col=c("blue", "orange"), pch=c(15,19), lty=c(3,1), lwd=c(2,2), bg="white")


colnames(chart_data) <- c("TPM", "Validation_Sample_Average", "Standard Deviation", sample_name)
write.table(chart_data, "chart.tsv", quote=FALSE, row.names=FALSE, sep="\t")

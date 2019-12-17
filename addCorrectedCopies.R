

args <- commandArgs(TRUE)
tempusFile <- args[1]
tumorPercent <- as.double(args[2])
geneFile <- args[3]

tempusData <- read.table(tempusFile, header = TRUE, sep="\t")

tempusData$TumorCorrectedCopies = ((2^(tempusData$Log2.Coverage.Ratio+1)-(2*(1-tumorPercent)))/tumorPercent)
#tempusData$TumorCorrectedCopies = (ifelse(tempusData$TumorCorrectedCopies > 0.1, tempusData$TumorCorrectedCopies, .1))
tempusData$RoundedCopies = round(tempusData$TumorCorrectedCopies,digits=0)

tempusData$Check = (tempusData$RoundedCopies==tempusData$Copy.Number)

geneData <- read.table(geneFile, header = TRUE, sep="\t")

first = TRUE
for( segment_row in 1:nrow(tempusData) ) {
	for( gene_row in 1:nrow(geneData) ) {
		if (geneData[gene_row, "Contig"] == tempusData[segment_row, "Chr"]) {
			if (geneData[gene_row, "END"] < tempusData[segment_row, "Start"]) {
				next
			} else if(tempusData[segment_row, "End"] < geneData[gene_row, "START"]) {
				next
			} else {
				if(first){
					tempusData[segment_row, "GENE"] = paste0(geneData[gene_row, "GeneSymbol"])
					first = FALSE
				} else if(is.na(tempusData[segment_row, "GENE"])){
					tempusData[segment_row, "GENE"] = toString(geneData[gene_row, "GeneSymbol"])
				} else {
					tempusData[segment_row, "GENE"] = paste0(geneData[gene_row, "GeneSymbol"], ", ", tempusData[segment_row, "GENE"])
				}
			}
		}
	}
}


write.table(tempusData, file="output.tsv", row.names=FALSE, sep="\t")

genes_only <- tempusData[complete.cases(tempusData[,"GENE"]),]

filtered_genes <- genes_only[genes_only$Copy.Number >= 5 | genes_only$Copy.Number == 0,]

write.table(filtered_genes, file="output2.tsv", row.names=FALSE, sep="\t")
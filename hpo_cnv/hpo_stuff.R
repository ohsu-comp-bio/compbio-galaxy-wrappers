
xlsxExists <- require("xlsx")

if (!xlsxExists){
	install.packages("xlsx")
	library("xlsx")
}

args <- commandArgs(TRUE)
case_dir <- args[1]
cnvGenesFile <- args[2]

hpo_files <- list.files(path=case_dir, pattern="*.xlsx")

first=TRUE
second=FALSE

for (file in hpo_files){
	ID <- unlist(strsplit(file, "\\."))[1]
	IDdata <- read.xlsx(paste0(case_dir, "/", file), sheetIndex=1)
	if (nrow(IDdata) < 1){
		next
	}
	#print(IDdata)
	IDdata <- as.character(IDdata[,"GENE_SYMBOL"])
	if (first){
		allCombined <- as.character(IDdata)
		names(IDdata) <- ID
		hpoData <- IDdata
		first=FALSE
		second=TRUE
		firstID <- ID
	} else {
		allCombined <- c(allCombined, IDdata)
		names(IDdata) <- ID
		hpoData <- cbind(IDdata, hpoData)
		colnames(hpoData)[1] <- c(ID)
		if (second){
			colnames(hpoData)[2] <- c(firstID)
			second=FALSE
		}
	}
}

hpoData <- apply(hpoData, 2, sort, decreasing=TRUE)

for (col in 1:(ncol(hpoData))){
	hpoData[duplicated(hpoData[,col]), col] <- NA
}

hpoData <- apply(hpoData, 2, sort, na.last=TRUE)

for (col in 1:(ncol(hpoData))){
	hpoData[is.na(hpoData[,col]), col] <- ""
}

write.table(hpoData, file=paste0(case_dir, "/HPOGeneLists.tsv"), row.names=FALSE, sep="\t")


deDuped <- unique(allCombined)

deDuped <- deDuped[order(deDuped)]

write.table(deDuped, file=paste0(case_dir, "/FullGeneList.tsv"), row.names=FALSE, sep="\t", quote=FALSE)

if (cnvGenesFile == "none"){
	print("no genes file found")
	cnvGenes <- as.data.frame(deDuped)
} else {
	cnvGenes <- read.table(cnvGenesFile, header=TRUE, sep="\t") 
	if (ncol(cnvGenes) > 1){
		print("Full genes tables found, extracting genes")
		cnvGenes <- as.data.frame(unlist(strsplit(as.character(cnvGenes[,c("GENE")]), ", ")))
		print(cnvGenes)
	}
}

colnames(cnvGenes) <- c("allCombined")

#write.table(allCombined, file="allCombined.tsv", sep="\t")
hpoCounts <- as.data.frame(table(allCombined))

hpoCounts <- hpoCounts[order(hpoCounts$Freq, decreasing=FALSE),]
#write.table(hpoCounts, file="hpoNums.tsv", row.names=FALSE, sep="\t")
hpo2 <- hpoCounts
hpo2 <- merge(hpoCounts, cnvGenes, by="allCombined" )

if (nrow(hpo2) < 1){
	print("No overlapping Genes")
}

hpo2 <- hpo2[order(hpo2$Freq, decreasing=TRUE),]

#write.table(hpo2, file="hpoNums2.tsv", row.names=FALSE, sep="\t")

#print(hpoData)

for (row in 1:nrow(hpo2)){
	freq = 0
	first=TRUE
	gene <- as.vector(hpo2[row, "allCombined"][[1]])
	#print(gene)
	for (col in 1:ncol(hpoData)){
		hpoID <- colnames(hpoData)[col]
		#print(hpoID)
		#print(gene)
		#print(as.vector(hpoData[,col]))
		#print(hpoData[,col])
		if (is.element(gene, as.vector(hpoData[,col]))){
			freq = freq + 1
			hpo2$Count[row] <- freq
			if (first){
				hpo2$HPO[row] <- paste0("HP:", hpoID)
				first = FALSE
			} else {
				hpo2$HPO[row] <- paste0(hpo2$HPO[row], ", HP:", hpoID)
			}
		}
	}
	
}
hpo2 <- hpo2[order(hpo2$Count, decreasing=TRUE),]

hpo2 <- hpo2[c("allCombined", "Count", "HPO")]

colnames(hpo2) <- c("Gene", "Num_HPO_Hits", "HPO_Terms")

write.xlsx(x=hpo2, file=paste0(case_dir, "/OverlappingGenesWithHPO.xlsx"), row.names=FALSE, sheetName=case_dir)

args <- commandArgs(TRUE)
cnvvcf <- args[1]
geneFile <- args[2]
readCountsFile <- args[3]
segmentsDir <- args[4]

library("tidyverse")
library("tidyr")
library("xlsx")

cnvData <- read.table(cnvvcf, header=TRUE, sep="\t", comment.char="", skip=11)
geneData <- read.table(geneFile, header=FALSE, sep="\t", skip=8)
readCounts <- read.table(readCountsFile, header=TRUE, sep="\t", comment.char="@")

colnames(geneData) <- c("ContigID", "ref", "type", "START", "END", "dot1", "sign", "dot2", "INFO")

#remove anything that isn't a gene or a region from the gene list
geneData <- geneData[geneData$type == "gene" | geneData$type == "region",]

#rename the sample column
colnames(cnvData)[10] <- "Nums"
#sperate the nums column
cnvData <- separate(cnvData, "Nums", c("Genotype", "CopyNumber", "Probes", "QCAll", "QCOne", "QCEnd", "QCStart"), sep=":", remove=FALSE, convert=TRUE)


#get rid of normals
cnvData <- cnvData[cnvData$Genotype != 0,]

#add an End column
cnvData$End <- gsub("END=", "", cnvData$INFO)
#make numeric
cnvData <- transform(cnvData, End = as.numeric(End), POS = as.numeric(POS))
geneData <- transform(geneData, END = as.numeric(END), START =  as.numeric(START))

#loop over each cnv
first=TRUE
contig = 0
for( segment_row in 1:nrow(cnvData) ) {
	#loop over the gene list
	for( gene_row in 1:nrow(geneData) ) {
		#if it is the first run through the gene list add the contig
		if (segment_row == 1 ) {
			#if the gene list is a region
			if ( geneData$type[[gene_row]] == "region" ){
				#get the info from the region 
				info <- as.vector(geneData$INFO[[gene_row]])
				infoSplit <- as.vector(strsplit(info, ";", fixed=TRUE))
				split <- as.vector(infoSplit[[1]])
				for( pair in split ) {
					pair <- toString(pair)
					#find the chromosome number and assign to contig
					if ( startsWith(pair, "chromosome=")){
						contig <- gsub("chromosome=", "", pair)
					}
				}
				#add the contig info to the CHROM column
				geneData[gene_row, "CHROM"] <- contig
			} else if ( geneData$type[[gene_row]] == "gene") {
				#if it is a gene add the contig info to the chrom column from the last region read
				geneData[gene_row, "CHROM"] <- contig
				#get the gene name
				info <- as.vector(geneData$INFO[[gene_row]])
				infoSplit <- as.vector(strsplit(info, ";", fixed=TRUE))
				split <- as.vector(infoSplit[[1]])
				for( pair in split ) {
					pair <- toString(pair)
					if ( startsWith(pair, "Name=")){
						#record the gene name in its own column
						geneData[gene_row, "GeneSymbol"] <- gsub("Name=", "", pair)
					}
				}
			}
		}
		#for every cnv row check if the chrom of the cnv and the chrom of the gene are the same, skip regions
		if ( geneData[gene_row, "CHROM"] == cnvData[segment_row, "X.CHROM"] & geneData$type[[gene_row]] == "gene" ) {
			#skip if the gene ends before the cnv starts
			if (geneData[gene_row, "END"] < cnvData[segment_row, "POS"]) {
				next
			#skip if the cnv ends before the gene begind
			} else if(cnvData[segment_row, "End"] < geneData[gene_row, "START"]) {
				next
			} else {
				#if this is the first gene cnv overlap create the GENE column
				if(first){
					cnvData[segment_row, "GENE"] = paste0(geneData[gene_row, "GeneSymbol"])
					first = FALSE
				#if it is the first gene overlap in this cnv segment, just add the gene
				} else if(is.na(cnvData[segment_row, "GENE"])){
					cnvData[segment_row, "GENE"] = toString(geneData[gene_row, "GeneSymbol"])
				} else {
					#if there is already a gene present, add it to the current list
					cnvData[segment_row, "GENE"] = paste0(geneData[gene_row, "GeneSymbol"], ", ", cnvData[segment_row, "GENE"])
				}
			}
		}
	}
	#for each cnv but not each gene
	#get the average count
	sum <- 0
	max <- 0
	min <- 100000
	inside <- FALSE
	for ( read_row in 1:nrow(readCounts)){
		if (readCounts$START[[read_row]] == cnvData$POS[[segment_row]] & readCounts$CONTIG[[read_row]] == cnvData$X.CHROM[[segment_row]]){
			inside <- TRUE
		}
		if (inside){
			readCount <- as.integer(readCounts$COUNT[[read_row]])
			sum <- sum + readCount
			max <- max(c(max, readCount))
			min <- min(c(min, readCount))
		}
		if (readCounts$END[[read_row]] == cnvData$End[[segment_row]]){
			inside <- FALSE
			break
		}
	}
	cnvData[segment_row, "AvgReadCount"] <- round(sum/as.integer(cnvData$Probes[[segment_row]]), digits=0)
	cnvData[segment_row, "MaxReadCount"] <- max
	cnvData[segment_row, "MinReadCount"] <- min
	print(paste0("Progress: ", cnvData[segment_row, "ID"], " : ", Sys.time()))
	
	#get the local freqs
	id <- paste0(cnvData[segment_row, "ID"])
	print(id)
	#Local freq same genotype
	genotype <- paste0(cnvData[segment_row, "Genotype"])
	grep_string <- paste0("grep -s ", id, " ", segmentsDir, "* | grep \t", genotype, ": | wc -l")
	localNum <- as.numeric(system(grep_string, intern=TRUE))
	print(localNum)
	numCase <- as.numeric(system(paste0("ls ", segmentsDir, " | wc -l"), intern=TRUE))
	print(numCase)
	localFreq <- as.numeric(localNum)/as.numeric(numCase)
	cnvData[segment_row, "LocalFreq"] <- round(localFreq, 3)
	
	#deletion frequency
	grep_string_del <- paste0("grep -s ", id, " ", segmentsDir,  "* | grep \t1: | wc -l")
	local_Del_Num <- as.numeric(system(grep_string_del, intern=TRUE))
	del_freq <- local_Del_Num/numCase
	cnvData[segment_row, "DelFreq"] <- round(del_freq, 3)
	
	#Dup freq
	grep_string_dup <- paste0("grep -s ", id, " ", segmentsDir, "* | grep \t2: | wc -l")
	local_Dup_Num <- as.numeric(system(grep_string_dup, intern=TRUE))
	dup_freq <- local_Dup_Num/numCase
	cnvData[segment_row, "DupFreq"] <- round(dup_freq, 3)
	print(localFreq)
}

#rename genotype del/dup
cnvData$Genotype <- ifelse(cnvData$Genotype == 1, "Del", "Dup")

#get size
cnvData$Size <- as.integer(cnvData$End) - as.integer(cnvData$POS)

#make pretty location column
cnvData$Location <- paste0("chr", cnvData$X.CHROM, ":", cnvData$POS, "-", cnvData$End)

#choose only relevant data
cnvData <- cnvData[c("Location", "GENE", "Genotype", "CopyNumber", "Probes", "AvgReadCount", "MaxReadCount", "MinReadCount", "Size", "LocalFreq", "DelFreq", "DupFreq", "QCAll", "QCOne", "QCStart", "QCEnd")]

#make a second pretty location
cnvData$Location2 <- gsub("chr", "", cnvData$Location)

#write to a tsv
#write.table(geneData, file="output_gene.tsv", row.names=FALSE, sep="\t")
write.table(cnvData, file="output_cnv.tsv", row.names=FALSE, sep="\t")

#write to an excel file
write.xlsx(x=cnvData, file="excel_cnv_output.xlsx", sheetName="AnnotatedCNVs", row.names=FALSE)

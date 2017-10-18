#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

##Check how many arguments were provided##
if (length(args) < 5) {
  stop("Not enough arguments supplied", call.=FALSE)
}

##Increase Java heap size##
options(java.parameters = "-Xmx6g")

##Load dependencies##
library(xlsx, quietly=TRUE)

##User-defined global variables##
mode = args[1]
file = args[2]
clinvar = args[3]
CSER_database = args[4]

if (mode == "cser") {
	ethnicity = "EUR"
	out = args[6]
	freq = 5

	keywords = c("coding-near-splice","codingComplex","synonymous-near-splice", "coding","frameshift", "frameshift-near-splice", "splice-donor", "splice-acceptor", "stop-gained", "stop-loss", "intron-near-splice", "missense-near-splice")

} else if (mode == "exome") {
	freq = args[5]
	gene_list = args[6]
	out = args[7]
	ethnicity = "Eur"

	keywords = c("coding-near-splice","missense","codingComplex","synonymous-near-splice", "coding","frameshift", "frameshift-near-splice", "splice-donor", "splice-acceptor", "stop-gained", "stop-loss", "intron-near-splice", "missense-near-splice")
}

###################################
#Function: main()
#Purpose: Run main set of commands
#Arguments: name of sample, ethnicity of sample, path to illumina carrierAll text file
###################################

main <- function () {
	##Import file##
	file <- read.csv(file, header=TRUE, sep="\t", stringsAsFactors=FALSE, na.strings=c("unknown","none",-100));

	if (mode == "exome") {
		gene_list <- read.csv(gene_list, header=FALSE, sep="\t");

		#Filter by gene list
		file[,"geneList"] <- sapply(strsplit(as.character(file$geneList), ","), "[", 1);
		merged <- merge(gene_list, file, by.y="geneList", by.x="V1");

		#Filter variants
		res <- filter_freq_refs(merged)
		res <- check_cser_database(res)
		res <- filter_quality(res)
		#res <- homozygotes(res)
		res <- match_clinvar(res)
	} else {
		#Filter variants
		res <- filter_freq_refs(file)
		res <- check_cser_database(res)
		res <- filter_quality(res)
		res <- homozygotes(res)
		res <- match_clinvar(res)
	}
	
	#Write output file
	write_excel(res)
}

###################################
#Function: filter_freq_refs
#Purpose: Filter variants based on frequency (5%) and available references
#Arguments: dataframe of illumina file
#Returns: resulting dataframe
###################################
filter_freq_refs <- function(data) {
	
	##filter variants on 1000Genomes and ExAC populations databases on a frequency determined by user
	if (mode == "exome") {
		data <- subset(data, data[paste("ExACExomes",ethnicity,".",sep="")] <= freq | is.na(data[paste("ExACExomes",ethnicity,".",sep="")]) == TRUE)
	}
	data <- subset(data, data[paste("X1000Genomes",ethnicity,".",sep="")] <= freq | is.na(data[paste("X1000Genomes",ethnicity,".",sep="")]) == TRUE)

	##filter variants on ClinVar or HGMD reference
	data <- subset(data, is.na(data$ClinVar) == FALSE | is.na(data$phenotypeHGMD) == FALSE | data$functionGVS %in% keywords)
	return(data)
}

###################################
#Function: check_cser_database
#Purpose: Check variants that have been previously reported in CSER database
#Arguments: resulting dataframe from filter_freq_refs()
#Returns: mylist$data, mylist$unreported, mylist$reported
###################################
check_cser_database <- function(data) {
	#read in CSER database
	unreport_database <- read.xlsx2(CSER_database, 2, as.data.frame=TRUE, header=TRUE)
	reported_database <- read.xlsx2(CSER_database, 3, as.data.frame=TRUE, header=TRUE)

	#match variants to CSER database of unreported variants
	unreported_variants <- subset(data, data$X.chrom %in% unreport_database$Chrom & data$positionHg19 %in% unreport_database$position)
	nomatch <- rbind(data,unreported_variants)
	nomatch <- nomatch[!(duplicated(nomatch) | duplicated(nomatch, fromLast = TRUE)), ] #Remove duplicates

	#match variants to CSER database of reported variants
	reported_variants <- subset(nomatch, nomatch$X.chrom %in% reported_database$Chrom & nomatch$positionHg19 %in% reported_database$position)
	nomatch <- rbind(nomatch,reported_variants)
	nomatch <- nomatch[!(duplicated(nomatch) | duplicated(nomatch, fromLast = TRUE)), ] #Remove duplicates
	
	mylist <- list("data" = nomatch, "unreported" = unreported_variants, "reported" = reported_variants)
	return(mylist)
}

###################################
#Function: filter_quality
#Purpose: Check quality filter and append onto resulting list
#Arguments: resulting list from check_cser_database()
#Returns: data$data, data$unreported, data$reported, data$filterFail, data$filterPass
###################################
filter_quality <- function(data) {	
	#Filter out variants that fail filter (QDFail)
	filterFail <- subset(data$data, data$data$filter != "PASS" & data$data$filter != "SnpCluster")
	data$filterFail <- filterFail
	#Filter variants that pass filter (PASS)
	filterPass <- subset(data$data, data$data$filter == "PASS" | data$data$filter == "SnpCluster" | is.na(data$data$phenotypeHGMD) == FALSE)	
	data$filterPass <- filterPass
	return(data)
}

###################################
#Function: homozygotes
#Purpose: Check for homozygotes and append onto resulting list
#Arguments: resulting list from filter_quality)
#Returns: data$data, data$unreported, data$reported, data$filterFail, data$filterPass, data$heterozygous, data$homozygous
###################################
homozygotes <- function(data) {
	##main homozyogous (change variable later)##
	names(data)
	depth <- data$filterPass[,9] #Careful, this column is hard-coded
	depth <- as.data.frame(matrix(unlist(strsplit(depth,"[(,)]")), ncol=3, byrow=TRUE)) #These lines will error if reads map to 3 places i.e. 34(0,14,21)
	heterozygous <- subset(data$filterPass, as.numeric(as.character(depth$V2)) > 1)
	homozygous <- subset(data$filterPass, depth$V2 == 0 | depth$V2 == 1)
	data$heterozygous <- heterozygous
	data$homozygous <- homozygous
	return(data)
}

###################################
#Function: match_clinvar
#Purpose: Check variants for ClinVar Significance, separates "Benign" from "Not benign", and appends onto resulting list
#Arguments: resulting list from homozygotes()
#Returns: data$data, data$unreported, data$reported, data$filterFail, data$filterPass, data$heterozygous, data$homozygous, data$benign, data$main
###################################
match_clinvar <- function(data) {
	#Read in ClinVar database and grab GRCH37(hg19) data
	clinvar <- read.csv(clinvar, sep="\t", header=TRUE)
	clinvar.ch37 <- subset(clinvar, clinvar$Assembly == "GRCh37")

	#Creates dataframe of Clinvar accessions and associated significance
	split_rcv <- strsplit(as.character(clinvar.ch37$RCVaccession), ";")
	rcv.clin.sig <- data.frame(ClinVar = unlist(split_rcv), ClinicalSignificance = rep(clinvar.ch37$ClinicalSignificance, sapply(split_rcv, length)));

	#Matches variants to clinvar significance using accessions
	data$main <- data$filterPass
	data$main$ClinVar <- sapply(strsplit(data$main$ClinVar, '\\.') ,'[', 1)
	data$main <- merge(data$main, rcv.clin.sig, by.x="ClinVar", by.y="ClinVar", all.x=TRUE)[,union(names(data$main), names(rcv.clin.sig))]

	#Creates two dataframes of benign variants and not benign variants based on clinvar significance
	data$benign <- data$main[grep("benign", data$main$ClinicalSignificance, ignore.case = TRUE),]
	data$main <- data$main[grep("benign", data$main$ClinicalSignificance, ignore.case = TRUE, invert=TRUE),] 

	return(data)
}

###################################
#Function: add_header
#Purpose: Creates a header in excel
#Arguments: wb=workbook, ws=worksheet, value="Name of header", startRow 
###################################
add_header <- function(wb, ws, value, startRow=NULL) {
	style <- CellStyle(wb) + Font(wb,color="red")
	if (is.null(startRow)) {
		rows <- getRows(ws)
		startRow=length(rows)+1
	}
	#create line break
	rows <- createRow(ws,rowIndex=startRow)
	sheetLineBreak <- createCell(rows, 1)
	setCellValue(sheetLineBreak[[1,1]], " ")

	#create Header
	rows <- getRows(ws)
	startRow=length(rows)+1
	rows <- createRow(ws, rowIndex=startRow)
	sheetTitle <- createCell(rows, 1)
	setCellValue(sheetTitle[[1,1]], value)
	setCellStyle(sheetTitle[[1,1]], style)

}

###################################
#Function: add_database
#Purpose: Creates a database in excel
#Arguments: data=dataframe, ws=worksheet, col.names, row.names, startRow 
###################################
add_dataframe <- function (data, ws, col.names, row.names, startRow=NULL) {
	if (is.null(startRow)) {
		rows <- getRows(ws)
		startRow=length(rows)+1
	}
	addDataFrame(data, ws, startRow, startColumn=1, col.names=col.names, row.names=row.names)
}

###################################
#Function: write_excel
#Purpose: Creates an excel spreadsheet of filtered variants and saves R data
#Arguments: res=resulting data from match_clinvar 
###################################
write_excel <- function(res) {
	wb <- createWorkbook(type="xlsx")
	ws <- createSheet(wb, sheetName = "Sheet1")

	#main variants to analyze
	add_dataframe(res$main, ws, col.names=TRUE, row.names=FALSE, startRow=1)

	#Write reported variants
	add_header(wb, ws, value="#Reported variants")
	add_dataframe(res$reported, ws, col.names=FALSE, row.names=FALSE)

	#Write placeholder for variants removed after analysis
	add_header(wb, ws, value="#Removed after analysis")

	#Write ClinVar Benign
	add_header(wb, ws, value="#Benign")
	add_dataframe(res$benign, ws, col.names=FALSE, row.names=FALSE)	

	#Write Unreported variants from CSER database
	add_header(wb, ws, value="#Variants seen before")
	add_dataframe(res$unreported, ws, col.names=FALSE, row.names=FALSE)

	#Write homozygous variants
	if (mode == "cser") {
		add_header(wb, ws, value="#Homozygous")
		add_dataframe(res$homozygous, ws, col.names=FALSE, row.names=FALSE)
	}

	#Write filter fail variants
	add_header(wb, ws, value="#FilterFail")
	add_dataframe(res$filterFail, ws, col.names=FALSE, row.names=FALSE)
	
	#save Workbook and save Rdata
	saveWorkbook(wb, out)
	#save(file=paste(sample,".Rdata",sep=""), res)
}

##Call main funciton
main();


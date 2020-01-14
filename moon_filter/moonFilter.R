
#get the moon file
args <- commandArgs(TRUE)
moonFile <- args[1]
#read the moon file into a dataframe
moonData <- read.table(file=moonFile, header = TRUE, sep="\t", quote="\"")
#make sure the names dont have spaces in them by changing all the spaces to periods
names(moonData) <- gsub(" ",".",names(moonData))

#COMPOUND HETEROZYGOTES
#check if mother and father vcfs have been provided
if (length(args) > 2) {
	#get the vcfs
	fatherFile <- args[2]
	motherFile <- args[3]
	#turn the vcfs into dataframes
	fatherData <- read.table(file=fatherFile, header=FALSE, sep="\t", comment.char="#")
	motherData <- read.table(file=motherFile, header=FALSE, sep="\t", comment.char="#")
	
	#create the ID column for all three dataframes with chromosome and position
	moonData$ID <- paste0(moonData$Chr., ":", moonData$Position)
	motherData$ID <- paste0(motherData$V1, ":", motherData$V2)
	fatherData$ID <- paste0(fatherData$V1, ":", fatherData$V2)

	#Add the mother and father columns to the moonData by checking if each ID is also in the parent vcf
	moonData$Mother <- is.element(moonData$ID, motherData$ID)
	moonData$Father <- is.element(moonData$ID, fatherData$ID)
	
	#get all the heterozygous variants
	hets <- moonData[moonData$Zygosity == "heterozygous", ]
	#order the variants by gene name
	hets <- hets[order(hets$Gene),]
	#filter
	hets <- hets[hets$Effect != "intron_variant",]
	hets <- hets[hets$Effect != "synonymous_variant",]
	#hets <- hets[!(grepl("prime_UTR_variant", hets$Effect)),]
	hets <- hets[hets$ClinVar != "Benign",]
	hets <- hets[hets$ClinVar != "Benign, Likely benign",]
	hets <- hets[hets$ClinVar != "Likely benign",]
	
	#remove all the deNovos
	hets <- hets[(hets$Mother == FALSE | hets$Father == FALSE),]
	#remove all the variants that are the only variant in that gene
	hets <- subset(hets,duplicated(hets$Gene) | duplicated(hets$Gene, fromLast=TRUE))
	
	#get lists of those from Mother and from Father by mergeing based in ID and throwing out any variant not in both dataframes
	fromMother <- merge(hets, motherData, by="ID", all=FALSE)
	fromFather <- merge(hets, fatherData, by="ID", all=FALSE)
	#keep only the variants in genes that appear in both parental lists
	compHets <- hets[is.element(hets$Gene, fromMother$Gene) & is.element(hets$Gene, fromFather$Gene),]
}

#De Novos

deNovos <- moonData[moonData$De.Novo == "true", ]

deNovos <- deNovos[deNovos$Depth >= 8, ]

deNovos <- deNovos[deNovos$Genotype.Quality >= 40, ]

deNovos <- deNovos[deNovos$Diploid.frequency <= .01, ]

deNovos <- deNovos[deNovos$gnomAD.frequency <= .01, ]

deNovos <- deNovos[deNovos$Effect != "intron_variant",]

deNovos <- deNovos[deNovos$Effect != "synonymous_variant",]

deNovos <- deNovos[!(grepl("non_coding", deNovos$Effect)),]

deNovos <- deNovos[!(grepl("prime_UTR_variant", deNovos$Effect)),]


#ClinVar
clinVars <- moonData[grepl("athogenic", moonData$ClinVar), ]


#Mito
mitos <- moonData[moonData$Chr. == "MT", ]

mitos <- mitos[mitos$Effect != "synonymous_variant",]

#Truncating
truncs <- moonData[grepl("stop_gained", moonData$Effect), ]
starts <- moonData[grepl("start_lost", moonData$Effect),]
stop_losses <- moonData[grepl("stop_lost", moonData$Effect),]
splices <- moonData[grepl("splice", moonData$Effect),]
truncs <- rbind(truncs, starts, stop_losses, splices)
truncs <- unique(truncs)
truncs <- truncs[truncs$Effect != "splice_region_variant&intron_variant",]

#put all the lists of wanted variants together
white_list <- rbind(deNovos, clinVars, mitos, truncs)
#check if we also want comphets and if so add
if (length(args) > 2){
	white_list <- rbind(white_list, compHets)
}

#get rid of duplicates 
filter_list <- unique(white_list)


#apply final filters
filter_list <- filter_list[filter_list$Disorder != "No associated disorder", ]

#write the output
write.table(filter_list, file="output.tsv", row.names=FALSE, sep="\t")

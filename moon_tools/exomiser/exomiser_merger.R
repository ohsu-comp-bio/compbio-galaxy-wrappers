### An R script to take in the exomiser output, merge them together and sort.

#get the command line inputs
args <- commandArgs(TRUE)


#loop over the input files
count = 1
for(i in args) {
	#if it's the first input file then just create a new data frame
	if ( count == 1 ) {
		exomiser_data <- read.table(args[count], header = TRUE, sep = '\t', comment.char = "")
		count = count + 1
	}	else {
		#if it is not the first input file then create a new data frame and concatenate with the old one
		new_data <- setNames(read.table(args[count], header = TRUE, sep = '\t', comment.char = ""), names(exomiser_data))
		exomiser_data <- rbind(exomiser_data, new_data)
		count = count + 1
	}
}


print(names(exomiser_data))
#sort the data by EXOMISER_GENE_COMBINED_SCORE
exomiser_data <- exomiser_data[order(-exomiser_data$EXOMISER_GENE_COMBINED_SCORE),]
#get just the top ten variants
top_ten <- exomiser_data[c(1:10),]

#write the merged and sorted dataframe to a csv

write.csv(exomiser_data, file = "merged_sorted.csv")

#write the top ten data to a csv

write.csv(top_ten, file = "top_ten.csv")
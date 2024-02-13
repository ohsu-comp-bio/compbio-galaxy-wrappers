# Script to create KMT2A-PTD illustration
# C. Frerich Aug 2023
# Input files: 1. GTF from Ensembl containing exon numbers and genomic coordinates for plotting 2. Chimeric junction parser output .txt
# Outputs single .pdf illustration

# libraries ---------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(data.table)
library(ggtranscript) #for plotting the reference genes
library(tidyr)
library(ggrepel)
library(grid)
library(gridExtra)
library(rtracklayer)

# Filepaths
args <- commandArgs(trailingOnly=TRUE)
gtf_path <- args[1]
ptd_path <- args[2]
output_path <- args[3]

# input files -----------------------------------------------------------
#GTF file downloaded from Ensembl archive: https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/
#gtf_path <- "/Users/chongbe/Downloads/CandaceTool/Homo_sapiens.GRCh37.75.gtf.gz" #hg19
gtf <- rtracklayer::import(gtf_path) #load GTF file into R - slow. Maybe more efficient to pull specific transcripts of interest

# Data wrangling ----------------------------------------------------------
gtf <- gtf %>% dplyr::as_tibble() %>%
  mutate(exon_number = as.integer(exon_number)) #convert to tibble format for usability, make sure exons are integers

#Automatically get sample ID, not sure how to automate this. Leave for Galaxy peeps
sample.id <- "23KD-287P0003"

#read in datafile
#ptd.data <- read.delim("/Users/chongbe/Downloads/CandaceTool/23KD-287P0003_Galaxy101.txt")
ptd.data <- read.delim(ptd_path)

# data check --------------------------------------------------------------
# remove PML-RARA from table
ptd.data <- ptd.data %>%
  filter(!Fusion == "PMLe6-RARAi2")

#check if data file is empty, if so abort.
if(nrow(ptd.data)==0){
  print("No KMT2A-PTD chimeric junctions were parsed, it's okay, this happens frequently")
  quit(status=5)
}

# Make PTD illustration ---------------------------------------------------
#reference coordinates
KMT2A_ref_exons <- gtf %>%
  dplyr::filter(
    !is.na(gene_name),
    type == "exon",
    ccds_id == "CCDS55791", #take only appropriate transcript ID
    exon_number<=20,) %>%
  mutate(exon_type = "reference", name = "KMT2A-PTD")


#what exons do we need to plot?
isoforms <- ptd.data %>%
  distinct(Exon1, Exon2) %>%
  mutate(UID = paste(Exon1, "-", Exon2))

# Get the unique exon combinations, ignoring minor variation in breakpoints.
# reformat coordinates to stack exons without introns
coord.to.plot = list()
for(i in 1:nrow(isoforms)){
  coord.to.plot[[i]] <- rbind(KMT2A_ref_exons %>%
                  mutate(order = case_when(exon_number<= isoforms[i,]$Exon1 ~ exon_number,
                                           exon_number>isoforms[i,]$Exon1 ~ exon_number + (isoforms[i,]$Exon1- isoforms[i,]$Exon2+1))),
                  KMT2A_ref_exons %>%
                  filter(exon_number %in% (isoforms[i,]$Exon2:isoforms[i,]$Exon1)) %>%
                  mutate(exon_type = "insertion",
                          order = (isoforms[i,]$Exon1+1):(isoforms[i,]$Exon1+(isoforms[i,]$Exon1- isoforms[i,]$Exon2+1))))

  #set order according to above
  coord.to.plot[[i]] <- coord.to.plot[[i]][order(coord.to.plot[[i]]$order),]

  #make new coordinates to plot
  coord.to.plot[[i]] <- coord.to.plot[[i]] %>%
    mutate(end2 = cumsum(width),
           start2 = end2 - (width-1),
           joinedExons = isoforms[i,]$UID,
           isoform = paste0("KMT2A-PTD", sep="\n", isoforms[i,]$UID, sep="\n", "CCDS55791.1"))

}
#combine listed dataframes created above into single data.frame for plotting
coord.to.plot <- rbindlist(coord.to.plot)

#Extract some summary information for plotting
ptd.data <- ptd.data %>%
  mutate(joinedExons = paste(Exon1, "-", Exon2))

ptd.summary <- ptd.data %>%
  group_by(joinedExons) %>%
  summarise(Nisoforms = dplyr::n(),
            Count = sum(Count),
            NormCount = sum(NormCount)) %>%
  mutate(isoform = paste0("KMT2A-PTD", sep="\n", joinedExons, sep="\n", "CCDS55791.1")) %>%
  arrange(Count)

#add X coord needed for plotting
test <- coord.to.plot %>%
  group_by(isoform) %>%
  filter(exon_type == "insertion") %>%
  filter(exon_number == min(exon_number))
ptd.summary$xCoord <- test$start2[match(ptd.summary$joinedExons, test$joinedExons)]


#now to plot
pdf(output_path, width = 14, height =5)  #Let galaxy peeps come up with appropriate, informative name for .pdf

ggplot() +
  geom_range(data = coord.to.plot, aes(xstart = start2,
                                xend = end2,
                                y = joinedExons,
                                fill = exon_type)) +
  scale_fill_manual(values = c("reference" = "white", "insertion" = "salmon"))+
  scale_y_discrete(limits = ptd.summary$joinedExons)+
  geom_text(data = coord.to.plot, aes( #add exon numbers
    x = start2 + (width / 2),
    y = joinedExons,
    label = exon_number),
    size = 2.5) +
  geom_text(data = ptd.summary, aes( #add summary labels
    x = xCoord,
    y = joinedExons,
    fontface = "bold",
    label = paste0("There are ", Nisoforms, " unique KMT2A (CCDS55791.1) isoforms/breakpoints joining exons ", joinedExons, " with a combined read count of ",
                   Count, " raw junction reads and ", NormCount, " normalized junction reads")),
    nudge_y = 0.4,
    size = 3) +
    labs(title = paste0(sample.id, ": Summary of exons involved in partial tandem duplication event."),
       subtitle = paste0("Exons downstream of 20 are not included in visualization")) +
  theme_classic()+
  theme( #clean up look
    legend.position = "bottom",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank()

)

#add table
#make table to add to plot
grid.newpage()
table_for_plot <- ptd.data %>%
  select(Fusion, Count, NormCount)

table_for_plot$Left.Partner <- paste0("KMT2A ", "Exon ", ptd.data$Exon1, sep="\n",
                                      ptd.data$Chrom1, ":", ptd.data$Coord1, sep="\n",
                                      "Transcript: ", ptd.data$CCDS1, sep="\n",
                                      "Exon junction type: ", ptd.data$RefStatus1
                                      )

table_for_plot$Right.Partner <- paste0("KMT2A ", "Exon ", ptd.data$Exon2, sep="\n",
                                       ptd.data$Chrom2, ":", ptd.data$Coord2, sep="\n",
                                       "Transcript: ", ptd.data$CCDS2, sep="\n",
                                       "Exon junction type: ", ptd.data$RefStatus2)


table_for_plot$Sequence <- ptd.data$CombinedSeq
table_for_plot <- arrange(table_for_plot, desc(Count) )


grid.table(table_for_plot, rows =NULL)

dev.off()

















































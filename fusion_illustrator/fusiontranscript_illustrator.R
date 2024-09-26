# Code to create fusion transcript illustrations
#C.Frerich 2024

# This code:
# Extracts coordinates for region of interest from ENSEMBL for both reference genes involved in the fusion event
# Extracts junction read coordinates from STAR-Fusion file
# Combines reference and fusion coordinates then rescales everything together so it will plot nicely
# Predicts protein coordinates included in the fusion from genomic coordinates
# Plots circos plot, fusion transcripts mapped to reference, predicted fusion protein, produces summary table
# # Code used to create input files May2024:
# # 1. Chimeric read information output from STAR-Fusion from two steps in galaxy
# # 2. V75 (hg19) ENSEMBL database created with the following code:

# Libraries ---------------------------------------------------------------
#data organization
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tidyr))
#for plotting
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtranscript)) #for plotting the reference genes
suppressPackageStartupMessages(library(ggbreak))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(ggbump))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(chimeraviz)) #for circos plot
suppressPackageStartupMessages(library(drawProteins)) #for plotting proteins
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(gridExtra))
#for mapping coordinates
suppressPackageStartupMessages(library(ensembldb)) #to get coordinates and stuff, this guy conflicts with dyplr a lot :(
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v75)) #this is hg19
suppressPackageStartupMessages(library(plyranges)) #need for protein info section
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(UniProt.ws)) #to get protein domains
suppressPackageStartupMessages(library(GenomicFeatures))
#misc
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(knitr))

# Load helper functions
source('ensembl_func.R')

options(datatable.rbindlist.check="warning")
options(datatable.optimize=1)
options(error=traceback)

# Input paths
args <- commandArgs(trailingOnly=TRUE)
sample.id <- args[1]
gtf_path <- args[2]
protein_ref <- args[3]
biomart_idmap <- args[4]
starfusion_results <- args[5]
starfusion_modified_results <- args[6]
# Output path
output <- args[7]


# LOAD FILES ----------------------------------------------------
#hg19 ensembl build, still need this to use genomeToProtein() which maps genomic coordinates to protein coordinates
edb <- EnsDb.Hsapiens.v75

#Ensembl database for transcript information
gtf <- import(gtf_path)
gtf <- gtf %>% #convert to Tibble for easy editing
  dplyr::as_tibble()

# #Ensembl database to look up length of protein for shading
# prts <- readRDS("D:/Heme_Fusion/Fusion illustration/ref_inputs/hg19_EnsDbHsapiensv75_proteins.RDS")
#Ensembl database to look up length of protein for shading
prts <- read.csv(protein_ref)

#Biomart download to look up Uniprot ID for specific ENSTID
uniprot_lookup <- read.csv(biomart_idmap, na.strings=c("","NA"))

#STAR-Fusion results for circos plot, has to only contain original star fusion columns or command will not recognize the file.
circos.fusions <- import_starfusion(starfusion_results,
                             "hg19")

#Import modified STAR-Fusion output file, with exons and other stuff
fusions <- read.delim(starfusion_modified_results,
                    header = TRUE,
                    check.names = TRUE,
                    sep = "\t")

# Fusion df cleaning ------------------------------------------------------
#remove some columns because they irritate me
fusions <- fusions %>%
  dplyr::select(-c(FUSION_CDS, PFAM_LEFT, PFAM_RIGHT))

#do some data wrangling for future steps
#Script requires the following columns:
#FusionName,
#Left.Gene, Left.Position,
#Right.Gene, Right.Position,
#Norm.JunctionRead, ENST_LEFT_ID, ENST_RIGHT_ID,
fusions <-  fusions %>%
  dplyr::rename(FusionName = name,
                Left.Gene = HGVSGene1,
                Left.Chr = chrom1,
                Left.Position = end1,
                Left.Exon = EXON_LEFT,
                Right.Gene = HGVSGene2,
                Right.Chr = chrom2,
                Right.Position = end2,
                Right.Exon = EXON_RIGHT,
                Norm.JunctionRead = NormalizedFrags) %>%
  dplyr::mutate(Left.Position = as.integer(Left.Position),
                Right.Position = as.integer(Right.Position),
                Left.Exon = as.integer(Left.Exon), #make sure these are numeric, will give warning, can ignore
                Right.Exon = as.integer(Right.Exon)) #make sure these are numeric



fusions <- split(fusions, f = fusions$FusionName)
#sort list so that fusion with the most (sum of all breakpoints) be on top
fusions <- fusions[order(sapply(fusions, function(i)colSums(i['Norm.JunctionRead'])), decreasing = TRUE)]

#determine the number of loops to perform
if(length(names(fusions)) >= 5){
  nloops = 5
} else {
  nloops = length(names(fusions))
}


# Reference genes ----------------------------------------------
# get gene coordinates from ENSEMBL to map the reference genes
ref_gene_coord <- list()
ref_gene_data <- list()
for(i in names(fusions[1:nloops])) {
  print(i)

  #pull out left and right gene names for the specific list you are on
  left_gene_name <- strsplit(i, "--")[[1]][1]
  right_gene_name <- strsplit(i, "--")[[1]][2]

  #pulls out transcript ID given by STAR-fusion, without the version info... because that is not included in ENSEMBL v75
  #note this takes the FIRST ENST ID if there are multiple
  left_ENST <-getLeftENST_function(fusions[[i]])
  right_ENST <- getRightENST_function(fusions[[i]])

  ref_gene_data[[i]] <- gtf %>%
    filter(gene_name %in% c(left_gene_name, right_gene_name)) %>%
    dplyr::mutate(exon_idx = as.integer(exon_idx)) #make sure these are numeric)

  #Get exon coordinates from reference genome for rescaling
  ref_gene_coord[[i]] <- ref_gene_data[[i]] %>%
    dplyr::filter(
      !is.na(gene_name),
      tx_id %in% c(left_ENST, right_ENST),
      case_when(#limit to exons of interest for each gene adding two upstream and downstream as buffer
                #Filter left gene exons
                gene_name == left_gene_name & nrow(fusions[[i]]) != sum(is.na(fusions[[i]]$Left.Exon))
                ~ exon_idx >= min(fusions[[i]]$Left.Exon, na.rm=T) - 2 &
                  exon_idx <= max(fusions[[i]]$Left.Exon, na.rm=T) + 2,

                #Filter right gene exons
                gene_name == right_gene_name & nrow(fusions[[i]]) != sum(is.na(fusions[[i]]$Right.Exon))
                ~ exon_idx >= min(fusions[[i]]$Right.Exon, na.rm=T) - 2 & #this seems to work even if it throws out a ridiculous negative number
                  exon_idx <= max(fusions[[i]]$Right.Exon, na.rm=T) + 2,

                #If there are only NA take all exons,NA will still throw warning, is okay to ignore
                TRUE
                ~ exon_idx >= min(ref_gene_data[[i]]$exon_idx, na.rm=T) &
                  exon_idx <= max(ref_gene_data[[i]]$exon_idx, na.rm=T)

                )) %>%
    dplyr::select(exon_seq_start,
                  exon_seq_end,
                  gene_name,
                  exon_idx) %>%
    dplyr::rename(start = exon_seq_start,
                  end = exon_seq_end,
                  exon_number = exon_idx)

   # # #do you wanna see them?
  # temp <- ref_gene_coord[[i]] %>% dplyr::filter(gene_name == left_gene_name)
  #
  # ggplot(data = temp, aes(xstart = start,
  #                                    xend = end,
  #                                    y = gene_name)) +
  #   geom_range(aes(fill = gene_name)) +
  #   geom_intron(data = to_intron(temp, "gene_name"))

  #melt reference genes for rescaling
  #ref_gene_coord[[i]] <- reshape2::melt(ref_gene_coord[[i]], id.var = c("gene_name", "exon_number"))

  ref_gene_coord[[i]] <- pivot_longer(ref_gene_coord[[i]], !c("gene_name", "exon_number"), names_to = "variable")
  ref_gene_coord[[i]] <- ref_gene_coord[[i]] %>%
    dplyr::mutate(value = as.integer(value)) #make sure these are numeric)


}


# Rescale -----------------------------------------------------------------
# Since both reference genes are on completely different chromosomes with wildly different coordinates
# we are going to rescale them (plus the fusion reads) to an arbitrary scale for plotting
temp.fusions <- list()
for(i in names(fusions[1:nloops])) {
  print(i)

  temp.fusions[[i]] <- bind_rows(
    dplyr::select(fusions[[i]], gene_name = Left.Gene, value = Left.Position),
    dplyr::select(fusions[[i]], gene_name = Right.Gene, value = Right.Position)
  ) %>%
    dplyr::mutate(variable = "start", exon_number = "fusion")

  pdata <- rbind(ref_gene_coord[[i]], temp.fusions[[i]]) #combine reference with experimental fusions
  #rescale everything together, they are added together in this fashion to ensure all parts scale correctly relative to each other
  pdata <-  pdata %>%
    group_by(gene_name) %>%
    mutate(new_value = rescale(value, to=c(0,100))) %>%
    ungroup()

  #extract ref gene coordinates and reformat for plotting
  ref_pdata <- pdata %>% dplyr::filter(!exon_number == "fusion")
  ref_gene_coord[[i]] <- pivot_wider(ref_pdata, id_cols = c("gene_name", "exon_number"), names_from = variable,
                                 values_from = new_value, values_fill = NA) #duplicates will throw error here
  ref_gene_coord[[i]]$strand <- ref_gene_data[[i]]$strand[match(ref_gene_coord[[i]]$gene_name, ref_gene_data[[i]]$gene_name)] #add strand information for plotting


  #extract junction information and and reformat for plotting
  junctions <- pdata %>% dplyr::filter(exon_number == "fusion")

  fusions[[i]] <- fusions[[i]] %>% #remap new scaled values from junctions to original fusions file for plotting
    dplyr::mutate(scaled.Left.Position = junctions$new_value[match(fusions[[i]]$Left.Position, junctions$value)],
           scaled.Right.Position = junctions$new_value[match(fusions[[i]]$Right.Position, junctions$value)],
           scaled.Junctions = case_when(fusions[[i]]$Norm.JunctionRead < 0 ~ 0.00, #just in case
                                        dplyr::between(fusions[[i]]$Norm.JunctionRead, 0, 10) ~ 0.01,
                                        dplyr::between(fusions[[i]]$Norm.JunctionRead, 11, 50) ~ 1,
                                        dplyr::between(fusions[[i]]$Norm.JunctionRead, 51, 100) ~ 3,
                                        fusions[[i]]$Norm.JunctionRead >= 100 ~ 9),
           y = 1.75, yend = 1.25)

  fusions[[i]] <- fusions[[i]] %>%
    dplyr::mutate(label.coord = abs(fusions[[i]]$scaled.Left.Position-fusions[[i]]$scaled.Right.Position)/2 + fusions[[i]]$scaled.Left.Position)

}

# Protein info ------------------------------------------------------------
# Pull information for protein plot
# Looks up the Uniprot ID for a specific ENSTID -- THIS a bit UNRELIABLE :(
uniprot_list <- list()
protplot_shading <- list()
for(i in names(fusions[1:nloops])) {
  print(i)

  #pull out left and right gene names for the specific list you are on
  left_gene_name <- strsplit(i, "--")[[1]][1]
  right_gene_name <- strsplit(i, "--")[[1]][2]

  #pulls out transcript ID given by STAR-fusion, without the version info... because that is not included in ENSEMBL v75
  left_ENST <-getLeftENST_function(fusions[[i]])
  right_ENST <- getRightENST_function(fusions[[i]])

  #Use ENST ID to look up the Uniprot ID using the uniprot lookup list
  uniprot_list[[i]] <- uniprot_lookup %>%
    dplyr::as_tibble() %>%
    filter(Transcript.stable.ID %in% c(left_ENST, right_ENST)) %>%
    distinct(Transcript.stable.ID, UniProtKB.Swiss.Prot.ID) %>%
    arrange(match(Transcript.stable.ID, c(left_ENST, right_ENST))) %>% #order table for easy graphing later
    pull(UniProtKB.Swiss.Prot.ID)

  # hack to bring in BCR protein for p190, otherwise end up with no plot.
  if(left_ENST == "ENST00000398512"){
    uniprot_list[[i]][2] <-  uniprot_list[[i]][1]
    #uniprot_list[[i]][1] = "A9UF05" #this is technically the matching uniprot ID, but star-fusion selected a short/dumb transcript
    uniprot_list[[i]][1] = "P11274" #this is the usual protein we want to plot
  }

  if(length(uniprot_list[[i]]) == 2 & any(is.na(uniprot_list[[i]])) == FALSE) {#if both uniprot ID were found, then find protein coordinates so we can add shading to the plot
    # get coordinates to shade portion of protein included in fusion
    #first use ENSP ID to find the end of the right protein for shading purposes
    right.ensp <- gtf %>% #get ENSP ID
      dplyr::filter(
        !is.na(gene_name),
        tx_id %in% c(right_ENST)) %>%
      distinct(protein_id) %>%
      pull(protein_id)

    right.ensp <-  prts %>% #get protein length
      dplyr::filter(names == right.ensp) #keep eye out for errors here


    #Now look up amino acid number using genomic coordinates
    #first pull out breakpoint coordinates, format them as granges object,
    # then use genometoprotein function which uses ENSEMBL hg19 to map
    prot_plot_shading <- fusions[[i]] %>%
      dplyr::select(Left.Chr, Left.Position, Right.Chr, Right.Position)

    prot_plot_shading <- prot_plot_shading %>% dplyr::select(Left.Chr, Left.Position) %>%
      bind_rows(prot_plot_shading[3:4] %>%
                  `colnames<-` (colnames(prot_plot_shading[1:2]))) %>%
      dplyr::rename(chrom = Left.Chr, start = Left.Position) %>%
      dplyr::mutate(end = start, chrom = gsub("chr", "", chrom))

    prot_plot_shading <- as(prot_plot_shading, "GRanges")

    #getting amnio acid number seems reliable
    prot_plot_shading <- genomeToProtein(prot_plot_shading, edb) #access Ensbl to map genomic coordinates to amino acid number
    prot_plot_shading <- as.data.frame(prot_plot_shading@unlistData)

    #now we have genomic coordinates mapped to transcript and amino acid
    prot_plot_shading <-  prot_plot_shading %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      dplyr::filter(tx_id %in% c(left_ENST, right_ENST))

    #Add the rest of the coordinates needed for plotting shading
    protplot_shading[[i]] <- prot_plot_shading %>%
      dplyr::mutate(ymin = case_when(tx_id == left_ENST ~ 1.6,
                                     tx_id == right_ENST ~ 0.6),
                    ymax = case_when(tx_id == left_ENST ~ 2.4,
                                     tx_id == right_ENST ~ 1.4),
                    start = case_when(tx_id == left_ENST ~ 0, TRUE ~ start),
                    end = case_when(tx_id == right_ENST ~ right.ensp$width, TRUE ~ end))
  }
}

# Plot --------------------------------------------------------------------
pdf(output, width = 12, height =8)

#Plot only junctions above the 2xSD background - note hard coded threshold
temp <- sapply(circos.fusions, function(x) x@split_reads_count >= 22.5)
circos.fusions <- circos.fusions[temp]

if(length(circos.fusions) >= 1) {
  plot_circle(circos.fusions) # make circos plot of fusions, works like a charm
  title(paste0(sample.id, " Circos plot: Limited to breakpoints with more than 22.5 raw junction reads"), line = 0) #add title
} else {
  ggplot() +
    annotate("text",
             x=1, y=1, size =4, colour = "red",
             label = paste0("No breakpoints in ", sample.id, " had more than 22.5 raw junction reads")) +
    theme_void()
}


for(i in names(fusions[1:nloops])) {

  print(i)

  #pull out gene name and transcript
  left_gene_name <- strsplit(i, "--")[[1]][1]
  right_gene_name <- strsplit(i, "--")[[1]][2]
  left_ENST <-getLeftENST_function(fusions[[i]])
  right_ENST <- getRightENST_function(fusions[[i]])

  #pull out total exon number for each transcript for labels
  left.tot.num.exons <- ref_gene_data[[i]] %>%
    dplyr::filter(
      !is.na(gene_name),
      tx_id %in% c(left_ENST)) %>%
    summarize(max(exon_idx))

  right.tot.num.exons <- ref_gene_data[[i]] %>%
    dplyr::filter(
      !is.na(gene_name),
      tx_id %in% c(right_ENST)) %>%
    summarize(max(exon_idx))

  #Save labels for transcript plot
  left_transcript_label <-  paste0(left_gene_name, sep="\n",
                                   left_ENST, sep="\n",
                                   fusions[[i]]$REFSEQ_LEFT_ID[[1]][1], sep="\n",
                                   fusions[[i]]$CCDS_LEFT_ID[[1]][1], sep="\n",
                                   "Total num exons:", left.tot.num.exons)

  right_transcript_label <-  paste0(right_gene_name, sep="\n",
                                    right_ENST, sep="\n",
                                    fusions[[i]]$REFSEQ_RIGHT_ID[[1]][1], sep="\n",
                                    fusions[[i]]$CCDS_RIGHT_ID[[1]][1], sep="\n",
                                    "Total num exons:", right.tot.num.exons)

  #make transcripts plot
  p = ggplot(data = ref_gene_coord[[i]], aes(xstart = start, xend = end, y = gene_name)) +
    geom_range(aes(fill = gene_name)) +
    geom_intron(data = to_intron(ref_gene_coord[[i]], "gene_name"), aes(strand = strand)) + #plot introns, with designated arrow direction
    #set X axis scale to plot left gene first
    scale_y_discrete(limits = c(right_gene_name, left_gene_name),
                     labels = c(right_transcript_label, left_transcript_label)) +
    #add exon numbers
    geom_text(data = dplyr::filter(ref_gene_coord[[i]], gene_name == left_gene_name),
              aes(x = start, y = gene_name, label = exon_number), nudge_y = 0.4) +
    geom_text(data = dplyr::filter(ref_gene_coord[[i]], gene_name == right_gene_name),
              aes(x = start, y = gene_name, label = exon_number), nudge_y = -0.4) +
    #add junctions
    ggbump::geom_sigmoid(
      data = fusions[[i]],
      aes(
        x = scaled.Left.Position,
        xend = scaled.Right.Position,
        y = y,
        yend = yend,
        group = rownames(fusions[[i]]),
        linewidth = scaled.Junctions),
      inherit.aes = F,
      direction = "y") +
    scale_linewidth_identity() + #prevent rescaling junction linewidth
    #Misc formatting
    labs(title = paste0(sample.id, ": normalized junction reads plotted to reference transcripts")) +
    theme( #clean up look
      legend.position = "none",
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )

  #Make protein plot
  if(length(uniprot_list[[i]]) == 2 & any(is.na(uniprot_list[[i]])) == FALSE) { #only plots if both uniprot ID are found
    prot_data <- drawProteins::get_features(paste(uniprot_list[[i]][1], uniprot_list[[i]][2])) #download protein information as json format
    #if the download above fails the loops gets truncated and does not finish running.
    prot_data <- drawProteins::feature_to_dataframe(prot_data) #convert into data.frame for ggplot

    #make sure order of proteins matches order of fusion
    prot_data <- prot_data %>%
      mutate(order = case_when(accession %in% uniprot_list[[i]][1] ~ 2,
                               accession %in% uniprot_list[[i]][2] ~ 1))

    #plot each element individually, then print
    prot_plot <- draw_canvas(prot_data)
    prot_plot <- draw_chains(prot_plot, prot_data)
    prot_plot <- draw_domains(prot_plot, prot_data, label_domains = FALSE)
    prot_plot <- draw_motif(prot_plot, prot_data)
    #p <- draw_phospho(p, prot_data, size = 8)
    # prot_plot <- draw_regions(prot_plot, prot_data)

    prot_plot <- prot_plot +
      geom_rect(data = protplot_shading[[i]], aes(xmin= start, xmax= end, ymin = ymin, ymax= ymax), alpha=0.2) +
      labs(title = "EXPERIMENTAL: Regions predicted to be included in fusion protein are shaded in gray")+
      scale_y_continuous(breaks = 1:2, labels = c(uniprot_list[[i]][2], uniprot_list[[i]][1]))+
      theme(legend.position="bottom", legend.text= element_text(size=8),
            axis.ticks = element_blank(), legend.title = element_blank(),
            axis.title.x = element_text(size=8), plot.title = element_text(size=8))

    } else {
    prot_plot <- ggplot() +
      annotate("text",
               x=1, y=1, size =4, colour = "red",
               label = "Error: could not successfully map both STAR-fusion selected transcripts to a Uniprot ID.
               \n Alternative protein coding transcripts with a Uniprot ID may exist") +
      theme_void()
  }

  #make table to add to plot
  table_for_plot <- fusions[[i]][ , c("FusionName", "PROT_FUSION_TYPE")]

  table_for_plot$Left.Partner <- paste0(fusions[[i]]$Left.Gene, sep="\n",
                                        fusions[[i]]$Left.Chr, ":", fusions[[i]]$Left.Position, sep="\n",
                                        "Transcript: ", fusions[[i]]$ENST_LEFT_ID, sep = ", ", fusions[[i]]$CCDS_LEFT_ID, sep = ", ", fusions[[i]]$REFSEQ_LEFT_ID, sep="\n",
                                        "Exon ", fusions[[i]]$Left.Exon)

  table_for_plot$Right.Partner <- paste0(fusions[[i]]$Right.Gene, sep="\n",
                                         fusions[[i]]$Right.Chr, ":", fusions[[i]]$Right.Position, sep="\n",
                                         "Transcript: ", fusions[[i]]$ENST_RIGHT_ID, sep = ", ", fusions[[i]]$CCDS_RIGHT_ID, sep = ", ", fusions[[i]]$REFSEQ_RIGHT_ID, sep="\n",
                                         "Exon ", fusions[[i]]$Right.Exon)


  table_for_plot$Junction.Reads <- fusions[[i]]$JunctionReadCount
  table_for_plot$Norm.Reads <- fusions[[i]]$Norm.JunctionRead
  table_for_plot <- head(table_for_plot, n = 3) #limit to only the top entries, otherwise will not fit on page

  table_for_plot <- tableGrob(table_for_plot, rows = NULL,
                              theme = ttheme_default(base_size = 7))

  #try this to add title
  grid.arrange(arrangeGrob(p),
               arrangeGrob(prot_plot),
               arrangeGrob(table_for_plot, top = 'Top 3 fusions'),
               nrow=3,
               heights = c(0.3, 0.4, 0.3),
               as.table=TRUE)

}

dev.off()

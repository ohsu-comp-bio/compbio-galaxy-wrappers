#########
# Generate abundance tables and metadata tables for downstream analysis
# 9-11-24
# Author: Christopher Suciu

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(singscore))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(cqn))
suppressPackageStartupMessages(library(scales))

# Input Filepaths
args <- commandArgs(trailingOnly=TRUE)
reference_path <- args[1]
print(reference_path)
counts_path <- args[2]
new_counts_path <- args[3]
meta_path <- args[4]
new_meta_path <- args[5]
# Output
output_uhr1 <- args[6]
output_uhr2 <- args[7]
output_tumor <- args[8]

# Load data
gene.info <- get(load(reference_path)) #gene.info
gene.info <- filter(as_tibble(gene.info), is.na(target_length)==F)
count.mat <- readr::read_delim(counts_path)
new_counts <- readr::read_delim(new_counts_path)
new_counts <- new_counts[,-1]
count.mat <- count.mat %>% left_join(., new_counts, by="hgnc_symbol")
count.mat <- count.mat %>% column_to_rownames("hgnc_symbol")
meta <- readr::read_delim(meta_path)
meta_new <- readr::read_delim(new_meta_path)
meta <- rbind(meta, meta_new)
meta <- mutate(meta, barcode=Run_Case_ID)

tumor.meta <- filter(meta, grepl("UHR", Cases)==F)
#We are also going to choose one case (earliest by RunID) if multiple versions are present
tumor.meta <- tumor.meta %>% arrange(Txome_RunID) %>% filter(!duplicated(Cases))
dge.target.only <- DGEList(
  count.mat[gene.info$gene_id,tumor.meta$barcode],
  remove.zeros = T)
#adjust gene.info again to remove zero count genes
gene.info <- filter(gene.info, gene_id %in% rownames(dge.target.only))
#get the stable genes for stingscore
#below is adapted from singscore::getStableGenes

tumor.meta <- filter(meta, grepl("UHR", Cases)==F)
tumor.meta <- tumor.meta %>% arrange(Txome_RunID) %>% filter(!duplicated(Cases))
dge.target.only <- DGEList(
  count.mat[gene.info$gene_id,tumor.meta$barcode],
  remove.zeros = T)
gene.info <- filter(gene.info, gene_id %in% rownames(dge.target.only))
filt.dge <- DGEList(count.mat[gene.info$gene_id,], remove.zeros = T)

logRPKM <- rpkm(filt.dge, gene.length=gene.info$target_length,normalized.lib.sizes = FALSE, log=T)
gene.info <- filter(as_tibble(gene.info), is.na(cds_length)==F)
logRPKM <- logRPKM[gene.info$gene_id, ]

meta <- mutate(
  meta,
  run_date=lubridate::ymd(sapply(strsplit(Run_Case_ID, "_"), "[[", 1))
)

uhr.meta <- filter(meta, grepl("UHR", Cases))

#For the runs with multiple UHRs can we choose one?
##first compute correlations
uhr.cors <- cor(logRPKM[,uhr.meta$barcode])
range(uhr.cors[lower.tri(uhr.cors)])#0.9238222 0.9960958

uhr.meta <- uhr.meta %>% filter(Run_Case_ID != "210901_VH00579_4_AAAGNTVHV/RDM-UHR-20210719-05-08-1.txt")
uhr.meta <- uhr.meta %>% filter(Run_Case_ID != "210901_VH00579_4_AAAGNTVHV/RDM-UHR-20210719-01-04-1.txt")
uhr.meta <- uhr.meta %>% filter(Run_Case_ID != "210901_VH00579_4_AAAGNTVHV/RDM-UHR-20210818-01-04-1.txt")

uhr.meta <- uhr.meta %>% filter(Run_Case_ID != "191023_NS500390_0386_AHFWWWBGXC/RDM-UHR-NEWFRAG-20191003-1.txt")
uhr.meta <- uhr.meta %>% filter(Run_Case_ID != "231030_VH00579_140_AAF5CTKM5_RDM-UHR-20231025-1.txt")
uhr.meta <- uhr.meta %>% filter(Run_Case_ID != "220228_VH00579_19_AAAH3GNHV/RDM-UHR-20220222-01-04-1.txt")

multi.uhr.count <- group_by(uhr.meta, Txome_RunID) %>% summarize(n()) %>% filter(`n()` > 1)
multi.uhr <- uhr.meta %>% filter(Txome_RunID %in% multi.uhr.count$Txome_RunID)

multi.cors <- bind_rows(lapply(group_split(multi.uhr, Txome_RunID), function(x){
  tmp.mat <- uhr.cors[x$barcode,x$barcode]
  tmp.mat[upper.tri(tmp.mat, diag=T)] <- NA
  as_tibble(tmp.mat, rownames="barcode") %>%
    pivot_longer(cols=-barcode) %>%
    filter(is.na(value)==F) %>%
    mutate(Txome_RunID=x$Txome_RunID[1])
}))

multi.cors <- filter(multi.cors, value >.95)
#are any runids missing?
setdiff(multi.uhr$Txome_RunID, multi.cors$Txome_RunID) %>%
  length %>%
  {stopifnot(. == 0)}
#so from the remaining can choose a single representative for each runid
multi.cors <- filter(multi.cors, !duplicated(Txome_RunID))
uhr.meta <- filter(uhr.meta, Txome_RunID %in% multi.uhr$Txome_RunID == F)
multi.uhr <- inner_join(multi.uhr, multi.cors, by=c("barcode", "Txome_RunID"))
uhr.meta <- bind_rows(
  uhr.meta,
  multi.uhr
)

stopifnot(uhr.meta %>% select(Txome_RunID) %>% n_distinct() == nrow(uhr.meta))

uhr.cors <- cor(logRPKM[,uhr.meta$barcode])
range(uhr.cors[lower.tri(uhr.cors)])#0.9452419 0.9954402

uhr.list <- group_split(uhr.meta, Group)
names(uhr.list) <- sapply(uhr.list, function(x) x$Group[1])

###########

logRPKM_UHR_Lot_2 <- logRPKM[, uhr.list$UHR_Lot_2$barcode]
logRPKM_UHR_Lot_1 <- logRPKM[, uhr.list$UHR_Lot_1$barcode]
logRPKM_tumor <- logRPKM[gene.info$gene_id,tumor.meta$barcode]

# Write to CSV
write.csv(logRPKM_UHR_Lot_2, output_uhr1)
write.csv(logRPKM_UHR_Lot_1, output_uhr2)
write.csv(logRPKM_tumor, output_tumor)
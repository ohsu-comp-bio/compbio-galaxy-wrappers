# absolute paths should always be used when inputing data
#setwd("/Volumes/Data/academic/svmlab/ml/cancer_projects/aml/OHSU/data/")
args <- commandArgs(trailingONLY = TRUE)
drug_screen <- read.csv("Processed_Drug_Data_beginning_to_10-4-2014_FINAL.csv", stringsAsFactors = FALSE, sep=",", check.names = FALSE)
AUC_All <- tapply(drug_screen[,"Inhibitor Interpreted Result: Area under the curve"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)
IC50_All <- tapply(drug_screen[,"Inhibitor Interpreted Result: IC50"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)

drug_screen <- drug_screen[which(drug_screen[,"Heme Malignancy: Diagnosis"] == "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"),]
AUC_AML <- tapply(drug_screen[,"Inhibitor Interpreted Result: Area under the curve"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)
IC50_AML <- tapply(drug_screen[,"Inhibitor Interpreted Result: IC50"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)

BeatAML_DD_All <- list()
BeatAML_DD_All$AUC <- AUC_All
BeatAML_DD_All$IC50 <- IC50_All
BeatAML_DD_All$AUC_AUC_Pearson <- cor(AUC_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_All$AUC_AUC_Spearman <- cor(AUC_All, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DD_All$IC50_IC50_Pearson <- cor(IC50_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_All$IC50_IC50_Spearman <- cor(IC50_All, method = "spearman", use = "pairwise.complete.obs")
save(BeatAML_DD_All, file = "BeatAML_DD_All.RData")

BeatAML_DD_AML <- list()
BeatAML_DD_AML$AUC <- AUC_AML
BeatAML_DD_AML$IC50 <- IC50_AML
BeatAML_DD_AML$AUC_AUC_Pearson <- cor(AUC_AML, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_AML$AUC_AUC_Spearman <- cor(AUC_AML, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DD_AML$IC50_IC50_Pearson <- cor(IC50_AML, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_AML$IC50_IC50_Spearman <- cor(IC50_AML, method = "spearman", use = "pairwise.complete.obs")
save(BeatAML_DD_AML, file = "BeatAML_DD_AML.RData")

RNASeq <- read.csv("dset3_OHSU_raw_expression_counts.csv", stringsAsFactors = FALSE, sep="\t", row.names = "gene", check.names = FALSE)
RNASeq <- log2(t(RNASeq) + 1)

common_patients_All <- intersect(rownames(IC50_All), rownames(RNASeq))
common_patients_AML <- intersect(rownames(IC50_AML), rownames(RNASeq))

gene_names <- c("CTLA-4", "Galectin-9", "PD-L1", "PD-L2", "PD1", "TIM-3", "VISTA", "LAG3", "TNFRSF4", "CD3G", "ARG1", "ARG2")
gene_ids <- c("CTLA4", "LGALS9", "CD274", "PDCD1LG2", "PDCD1", "HAVCR2", "C10orf54", "LAG3", "TNFRSF4", "CD3G", "ARG1", "ARG2")
RNASeq_All <- RNASeq[common_patients_All, gene_ids]
RNASeq_AML <- RNASeq[common_patients_AML, gene_ids]
colnames(RNASeq_All) <- gene_names
colnames(RNASeq_AML) <- gene_names

AUC_All <- AUC_All[common_patients_All,]
IC50_All <- IC50_All[common_patients_All,]
#selected_drugs <- which(colSums(!is.na(AUC_All)) >= 6)
#AUC_All <- AUC_All[,selected_drugs]
#IC50_All <- IC50_All[,selected_drugs]
AUC_AML <- AUC_AML[common_patients_AML,]
IC50_AML <- IC50_AML[common_patients_AML,]
#selected_drugs <- which(colSums(!is.na(AUC_AML)) >= 6)
#AUC_AML <- AUC_AML[,selected_drugs]
#IC50_AML <- IC50_AML[,selected_drugs]

BeatAML_DG_All <- list()
BeatAML_DG_All$AUC <- AUC_All
BeatAML_DG_All$IC50 <- IC50_All
BeatAML_DG_All$RNASeq <- RNASeq_All
BeatAML_DG_All$AUC_RNASeq_Pearson <- cor(AUC_All, RNASeq_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_All$AUC_RNASeq_Spearman <- cor(AUC_All, RNASeq_All, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DG_All$IC50_RNASeq_Pearson <- cor(IC50_All, RNASeq_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_All$IC50_RNASeq_Spearman <- cor(IC50_All, RNASeq_All, method = "spearman", use = "pairwise.complete.obs")
save(BeatAML_DG_All, file = "BeatAML_DG_All.RData")

BeatAML_DG_AML <- list()
BeatAML_DG_AML$AUC <- AUC_AML
BeatAML_DG_AML$IC50 <- IC50_AML
BeatAML_DG_AML$RNASeq <- RNASeq_AML
BeatAML_DG_AML$AUC_RNASeq_Pearson <- cor(AUC_AML, RNASeq_AML, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_AML$AUC_RNASeq_Spearman <- cor(AUC_AML, RNASeq_AML, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DG_AML$IC50_RNASeq_Pearson <- cor(IC50_AML, RNASeq_AML, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_AML$IC50_RNASeq_Spearman <- cor(IC50_AML, RNASeq_AML, method = "spearman", use = "pairwise.complete.obs")
save(BeatAML_DG_AML, file = "BeatAML_DG_AML.RData")

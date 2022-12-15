    # VERSION: 0.9.4

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gtable))

source('dsp_inputs.R')
source('evaluation.R')
source('helpers.R')
source('normalization.R')
source('plots.R')
source('summarization.R')

options(datatable.rbindlist.check="warning")
options(datatable.optimize=1)
options(error=traceback)

### START ARGS
args <- commandArgs(trailingOnly=TRUE)
my_samp <- args[1]
runid <- args[2]
coh <- args[3]

# Metadata filepaths
datadir <- args[4]
ab_info <- args[5]
low_probes <- args[6]
control_type <- args[7]
dsp_meta <- args[8]
igg_info <- args[9]

# Set constants
exp.regex <- "[0-9]{8}"
clia_abs <- c("c-Myc", "Ki-67", "Cyclin E1", "Cyclin D1", "Cyclin B1", "ATR_pS428", "ATM (phospho S1981)",
              "PARP", "Cleaved Caspase 9", "CD95/Fas", "BIM", "BCLXL", "BCL6", "Bcl-2", "BAD",
              "Pan-AKT", "PTEN", "Phospho-PRAS40 (T246)", "PLCG1", "INPP4B", "Phospho-GSK3B (S9)", "Phospho-AKT1 (S473)",
              "pan-RAS", "p44/42 MAPK ERK1/2", "Phospho-p90 RSK  (T359/S363)", "MET", "Phospho-MEK1 (S217/S221)", "Phospho-p44/42 MAPK ERK1/2 (T202/Y204)", "HER2_p1248", "Her2", "EGFR", "Phospho-c-RAF (S338)", "BRAF",
              "p53", "Phospho-p38 MAPK (T180/Y182)", "PR", "Phospho-JNK (T183/Y185)", "ER-alpha", "Androgen Receptor", "ARID1A",
              "PD-L1", "CD8", "CD68", "CD4", "CD3", "CD20", "Beta-2-microglobulin")

# Load metadata
paths <- data.table(read.xlsx(ab_info, sheet="parsed"))
igg.info <- data.table(openxlsx::read.xlsx(igg_info))
exp.low <- data.table(openxlsx::read.xlsx(low_probes, colNames=F))
control.type <- data.table(openxlsx::read.xlsx(control_type, startRow=2))
stopifnot(control.type[,.N,by=.(lower_secondary, name, type)][,all(N==1)])
dsp.meta <- data.table(read.xlsx(dsp_meta))

# Check to see if my_samp is in dsp.meta, stop if it's not.
stopifnot(nrow(dsp.meta[`Specimen.ID` %in% my_samp]) > 0)

batch.dt <- data.table(files=list.files(datadir, pattern="[-0-9A-Za-z_ ]*nit[-0-9A-Za-z_ ]*ataset[-0-9A-Za-z_ ]*[0-9]{8}.xlsx", recursive = T, full.names=T))
batch.dt[,batch:=str_extract(files, exp.regex)]

res.list <- process_batches(batch.dt, sheet='Exported dataset')
qc.meta <- res.list$meta

#clean up sample labeling if necessary
qc.meta[`Segment (Name/ Label)` == "Segment 1",.(.N, max_roi=length(unique(croi))),by=.(sample_id, batch)][(N != 3) | (N != max_roi)]
qc.meta[sample_id=="01142", `:=`(sample_id="2001142")]
#Pull out the QC corrected data

all.abund <- Reduce(function(x,y){

  merge(x,y, by=c("ProbeName"))

}, lapply(res.list$data, function(x) x$qc[,-c(1:3)]))
stopifnot(nrow(all.abund) == res.list$data[[1]]$qc[,.N])
stopifnot((ncol(all.abund)-1) == sum(sapply(res.list$data, function(x) ncol(x$qc)-4 )))
abund.mat <- log2(as.matrix(all.abund[,-1,with=F]))
rownames(abund.mat) <- all.abund[,ProbeName]

# Get TMA submatrix
tma.meta <- qc.meta[`Segment (Name/ Label)`=="Full ROI" | `Segment (Name/ Label)`=="Geometric Segment"]
tma.meta[,lower_sample:=tolower(sample_id)]

stopifnot(control.type[,.N,by=.(lower_secondary, name, type)][,all(N==1)])

tma.meta <- merge(tma.meta, control.type[,.(lower_sample=lower_secondary, name, type)], by="lower_sample", all=F)
stopifnot(tma.meta[,.N,by=batch][,all(N==14)])
tma.abund <- abund.mat[,tma.meta$barcode]

## QC of experimental samples with respect to ROI
exp.meta <- qc.meta[`sample_id` %in% dsp.meta$Specimen.ID]
exp.meta$`Segment (Named/ Label)` <- ifelse(exp.meta$`Segment (Name/ Label)` == 'Geometric Segment', "Segment 1", exp.meta$`Segment (Name/ Label)`)
exp.meta$`Segment (Name/ Label)` <- ifelse(exp.meta$`Segment (Name/ Label)` == 'Full ROI', "Segment 1", exp.meta$`Segment (Name/ Label)`)

exp.abund <- abund.mat[, exp.meta$barcode]
#note we do this with the raw data as otherwise will constantly have to monitor due to normalization fluxes
exp.meta <- flag_roi(exp.abund, exp.meta, roi.thresh=.85)
#note the handful that didn't pass this QC, note that 432 would likely be considered marginal and kept
exp.meta[rm_croi == T,.(`Segment (Name/ Label)`, batch, sample_id, croi, max_cor, rm_croi)]
#for now, removing all
exp.meta <- exp.meta[rm_croi == F]
exp.abund <- abund.mat[, exp.meta$barcode]
#at this point save the metadata and raw abundance
save(exp.meta, exp.abund, tma.abund, tma.meta, file="tmp_dsp_batch1_11_expt_abund.RData")

### Normalization and summarization of cohort
# Determine relevant samples / batches and compute normalization factors
relevant.meta <- exp.meta[sample_id %in% dsp.meta$Specimen.ID]
# Figure out how we want to number batches.
relevant.meta[,num_batch:=batch]
relevant.abund <- exp.abund[,relevant.meta$barcode]
tma.meta <- tma.meta[batch %in% relevant.meta$batch]
# Figure out how we want to number batches.
tma.meta[,num_batch:=batch]
tma.abund <- abund.mat[,tma.meta$barcode]
cntrl.abs <- setdiff(rownames(tma.abund), exp.low[[1]])
segment.proc <- preprocess_dsp_tma(tma.meta, tma.abund, relevant.meta, relevant.abund, igg.map=igg.info, bg.method=c('none'), controls=cntrl.abs, use.type='quant', k=2, num.roi.avg=1)
#can save here

# Use the summarized metadata to deal with replicates and form groups
#create per segment summarized metadata
avg.meta <- rbindlist(lapply(segment.proc, function(x){
  x$meta
}), idcol="Segment (Name/ Label)")
## Batch 19 samples lost here
my.meta <- merge(dsp.meta[,.(num_batch=Date_run, sample_id=Specimen.ID, Specimen.ID, cohort, code)],
                 avg.meta[,.(`Segment (Name/ Label)`, num_batch, sample_id, avg_barcode)], by=c("num_batch", "sample_id"))
#minus two from above

#choose a technical replicate
## First determine how correlated they are to each other in segment abundance
sapply(names(segment.proc), function(x){

  tmp.meta <- my.meta[`Segment (Name/ Label)`==x]
  tmp.abund <- segment.proc[[x]]$avg_abund[,tmp.meta$avg_barcode]
  pair.cors <- sapply(split(tmp.meta, by=c("Specimen.ID", "code")), function(y){
    min(cor(tmp.abund[,y$avg_barcode]))
  })
  pair.cors
})

my.meta <- my.meta[!duplicated(cbind(`Segment (Name/ Label)`, Specimen.ID, num_batch, code)),]
#Here define reference vs experimental
my.meta[cohort==coh,Best_Response:="Ref"]
save(my.meta, segment.proc, file="tmp_proc_data.RData")

# Form antibody scores as the quantiles relevant to reference cohort
ref_samps <- my.meta[Best_Response == "Ref",unique(sample_id)]
ref_samps <- ref_samps[!ref_samps %in% my_samp]
quant.list <- score_abs(segment.proc, ref.samples=ref_samps ,score.type="quant")
#combine
pat.quants <- rbindlist(lapply(quant.list, "[[", "scores"), idcol="Segment (Name/ Label)")
my.scores <- merge(my.meta, pat.quants, by=c("Segment (Name/ Label)", "avg_barcode"))
# JHL: In the case of identical sample id's, get rid of the one that we are not currently analyzing.
my.scores <- my.scores[!(`sample_id` == my_samp & `num_batch` != runid)]

#getting pathways in order
use.paths <- paths[analysis_pathway %in% c("Expression Controls", "N/A")==F]
use.paths[analysis_pathway == "Tumor Markers", analysis_pathway:="Other Markers"]
use.paths <- use.paths[order(ab)]
path.ord <- c("Cell Cycle", "Cell Death", "PI3K/AKT pathway", "RAS/MAPK pathway", "Other Markers", "IO CLIA", "IO RUO")
use.paths[,`:=`(path_ord=factor(analysis_pathway, levels=path.ord, ordered=T),
                ab_ord=factor(ab, levels=ab, ordered=T))]

my.scores <- merge(use.paths[,.(ab_ord, ProbeName, path_ord)], my.scores, by="ProbeName", all=F)
my.scores[,comb_id:=Specimen.ID]
my.scores[,segment_label:=ifelse(`Segment (Name/ Label)` == "Segment 1", "tumor", "stroma")]
my.scores[,sample_ord:=sample_id]

pt.ord <- my.scores[,.N,by=.(comb_id, Best_Response)][order(Best_Response)]
my.scores[,patient_ord:=factor(comb_id, levels=pt.ord$comb_id, ordered=T)]
#output abundance for reference samples
ref.abund <- rbindlist(lapply(quant.list, "[[", "ref_abund"), idcol="Segment (Name/ Label)")
ref.abund[,segment_label:=ifelse(`Segment (Name/ Label)` == "Segment 1", "tumor", "stroma")]
ref.abund <- merge(use.paths[,.(ab_ord, ProbeName, path_ord)], ref.abund, by="ProbeName", all=F)
save(my.scores, ref.abund, file="tmp_results.RData")

# Steps used to provide data for antibody plots.
melt.tma <- data.table(reshape2::melt(tma.abund, as.is=T))
names(melt.tma) <- c("ProbeName", "barcode", "abundance")
melt.tma <- merge(melt.tma, tma.meta[,.(barcode, name, batch)], by="barcode")
stopifnot(length(setdiff(melt.tma$ProbeName, igg.info$ProbeName)) == 0)
melt.tma <- merge(igg.info[,.(ProbeName, igg)], melt.tma, by="ProbeName", all=T)
melt.tma[,fac_batch:=factor(batch)]
#add in values for corresponding igg
melt.tma[ProbeName %in% c("Ms IgG1",  "Ms IgG2a", "Rb IgG"), igg:=ProbeName]
melt.tma <- merge(melt.tma, melt.tma[ProbeName %in% c("Ms IgG1",  "Ms IgG2a", "Rb IgG"),.(igg=ProbeName, barcode, igg_abund=abundance)], by=c("igg", "barcode"), all.x=T, all.y=F)
melt.tma[,perc_igg:=(abundance/igg_abund)*100]
melt.tma[,ceil_igg:=pmin(perc_igg, 100)]

# Write out csv for Westgard rules script in Galaxy wf
write.csv(melt.tma, file=paste0(args[10]), row.names=F)

ref.batches <- melt.tma[batch != runid]
cur.batch <- melt.tma[batch == runid]
# Get the number of runs from the metadata sheet
run_no <- length(unique(dsp.meta$Date_run))
use.pal <- scales::hue_pal()(run_no)

# Remove unused factors
samp.scores <- my.scores[`sample_id` == my_samp]
samp.scores$patient_ord <- droplevels(samp.scores$patient_ord, except=my_samp)

# Overall Plots
plot.list <- loli_plot(score.dt=samp.scores, ref.dt=ref.abund)
pdf(file=paste0(args[11]), width=16, height=16)

for (tums in names(plot.list)){
  show(plot.list[[tums]])
}

# Create and write table of normalized counts.
tt <- ttheme_default(base_size = 16)
score_out <- samp.scores %>% select(ProbeName,segment_label,norm)
score_tum <- score_out[`segment_label` == 'tumor']
score_str <- score_out[`segment_label` == 'stroma']

g <- tableGrob(score_tum[1:38,1:3], rows = NULL, theme = tt)
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))

g1 <- tableGrob(score_tum[39:76,1:3], rows = NULL, theme = tt)
g1 <- gtable_add_grob(g1,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 2, b = nrow(g1), l = 1, r = ncol(g1))
g1 <- gtable_add_grob(g1,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, l = 1, r = ncol(g1))

haligned <- gtable_combine(g,g1, along=1)
grid.newpage()
grid.draw(haligned)

if (nrow(score_str) > 0) {
  g <- tableGrob(score_str[1:38,1:3], rows = NULL, theme = tt)
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, l = 1, r = ncol(g))

  g1 <- tableGrob(score_str[39:76,1:3], rows = NULL, theme = tt)
  g1 <- gtable_add_grob(g1,
                        grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 2, b = nrow(g1), l = 1, r = ncol(g1))
  g1 <- gtable_add_grob(g1,
                        grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 1, l = 1, r = ncol(g1))

  haligned <- gtable_combine(g,g1, along=1)
  grid.newpage()
  grid.draw(haligned)

}

for (ab in clia_abs){

  print(ab)
  q.plot <- ggplot(data=ref.batches[ProbeName == ab], mapping=aes(x=name, y=abundance)) +
    geom_boxplot(outlier.shape=NA) + geom_jitter(mapping=aes(color=fac_batch),linewidth=3, height=0, width=.15) +
    geom_jitter(data=cur.batch[ProbeName == ab], linewidth=3, width=.15, height=0) +
    scale_color_manual(values=setNames(use.pal, unique(dsp.meta$batch)), name="Previous Batches") +
    theme_bw() + xlab("") + ylab("log2 Abundance") + ggtitle(paste("Antibody: ", ab, "                Run ID: ", runid))

  plot(q.plot)

}

dev.off()
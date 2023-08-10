# VERSION: 1.2.1
# Version history
# 1.0.1 - Named all arguments, reformatting, allow for reference to contain only segment 1
# 1.0.2 - Use after_stat instead of stat(quantile) per deprecation msg, exclude comp samp from ref_samps
# 1.0.3 - Added tabular representation (.csv) of pairwise deltas
# 1.1.0 - Consolidated plots to one PDF file
# 1.2.0 - Removed selected_pos and pos_cntrls files, read all antibodies from GeoMx input to ensure we are using
#       - equivalent naming.  Split plots in to multiple pages for clarity and plot all antibodies.
# 1.2.1 - removed delta ridge plots and added cover sheet

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(gtable))

source('dsp_inputs.R')
source('evaluation.R')
source('helpers.R')
source('normalization.R')
source('plots.R')
source('summarization.R')

### START ARGS
args <- commandArgs(trailingOnly=TRUE)
my_samp <- args[1]
my_samp_compare <- args[2]
runid <- args[3]
runid_compare <- args[4]

# Metadata filepaths
ab_info <- args[5]

# dsp_runner intermediates
exp.out <- args[6]
seg.proc.out <- args[7]

# Output files
report.out <- args[8]

paths <- data.table(read.xlsx(ab_info, sheet="parsed"))
load(exp.out)
load(seg.proc.out)

## Pairwise Delta Distribution from RUV Normalized TMA
tma.meta[,tmpval:=1]
cl.repmat <- reshape2::acast(barcode~name, value.var="tmpval", data=tma.meta[type=="quant"], fun.aggregate = function(x) as.integer(length(x) > 0))
#and then to line 202 to create cl.pcs
cl.pcs <- compute_factors(tma.abund, cl.repmat)

use.abs <- paths[Manuscript=='y']

#Apply RUV
all.norm <- apply_norm(tma.abund, cl.pcs$pcs, k=2, controls=rownames(tma.abund))

#This is used to create comb.deltas and probe.mads
melt.dsp <- data.table(reshape2::melt(data=all.norm, as.is=T))
names(melt.dsp) <- c("ProbeName", "barcode", "value")
#Dan's code required date column which is used by altering line 161 of his code for any relevant data table that needs dates
tma.meta[,date:=lubridate::mdy(sapply(strsplit(Scan_ID, "_|\\s"), function(x) x[length(x)]))]
batch.ord <- tma.meta[,.N,by=.(date)]
batch.ord <- batch.ord[order(date)]
tma.meta[,date_fac:=factor(as.character(date), levels=as.character(batch.ord$date), ordered=T)]
melt.dsp <- merge(melt.dsp, tma.meta, by=c("barcode"))

by.ab <- split(melt.dsp, by="ProbeName")
val.range <- range(melt.dsp$value)

# THIS TABLE SHOWS STATS FOR NORMALIZED TMA's
cv.sum <- melt.dsp[,.(sd=sd(value), mean=mean(value), mad=mad(value),.N),by=.(ProbeName, name)]

#looks like stats are calculated for all TMA-ab combos
#But because more than one TMA per ab, MAD value for a given ab is taken from the worst of its associated tmas!!!!
# where is max.mad used, I don't see it anywhere
max.mad <- cv.sum[,.(max_mad=max(mad)),by=ProbeName][order(-max_mad)]

cv.sum[,N:=NULL]
names(cv.sum)[2] <- "cell_line"

delta.ab.pos <- rbindlist(lapply(split(melt.dsp, by=c("ProbeName", "name")), function(x){
  tmp.samps <- data.table(expand.grid(list(sample1=x$barcode, sample2=x$barcode), stringsAsFactors = F))[sample1 != sample2]
  tmp.deltas <- merge(tmp.samps, melt.dsp[ProbeName == x$ProbeName[1],.(ProbeName, sample1=barcode, sample1_value=value)], by="sample1", allow.cartesian=T)
  tmp.deltas <- merge(tmp.deltas, melt.dsp[ProbeName == x$ProbeName[1],.(ProbeName, sample2=barcode, sample2_value=value)], by=c("sample2", "ProbeName"))
  #removes redundant pairs, only keep x-->y not y-->x
  tmp.deltas[,idx:=sapply(strsplit(paste(sample1, sample2), " "), function(x) paste(sort(x), collapse=" "))]
  tmp.deltas <- tmp.deltas[!duplicated(idx)]
  tmp.deltas
}))
delta.ab.pos[,delta:=sample1_value-sample2_value]
comb.deltas <- rbind(
  cbind(delta.ab.pos, type="Positive Cell Lines")
)
comb.deltas[,ab_fac:=factor(ProbeName, levels=rev(paths$ProbeName), ordered=T)]

# Create pdf file to write to: PLOTTING STARTS HERE
pdf(file=paste0(report.out), width=16, height=16)

# Form antibody scores as the quantiles relevant to reference cohort
ref_samps <- my.meta[Best_Response == "Ref",unique(sample_id)]
ref_samps <- ref_samps[!ref_samps %in% my_samp]
ref_samps <- ref_samps[!ref_samps %in% my_samp_compare]
# If there is no ref data for segment 3 (sarcomas) then we'll just look at segment 1.
if (all(unique(my.meta[`sample_id` %in% ref_samps]$`Segment (Name/ Label)`) == "Segment 1")) {
  quant.list <- score_abs(segment.proc$`Segment 1`, ref.samples=ref_samps, stroma=F,score.type="quant")
  pat.quants <- quant.list$scores[,"Segment (Name/ Label)":="Segment 1"]
}else{
  quant.list <- score_abs(segment.proc, ref.samples=ref_samps, stroma=T,score.type="quant")
  #combine
  pat.quants <- rbindlist(lapply(quant.list, "[[", "scores"), idcol="Segment (Name/ Label)")
}
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
#at this point add in segment labels
my.scores[,segment_label:=ifelse(`Segment (Name/ Label)` == "Segment 1", "tumor", "stroma")]
my.scores[,sample_ord:=sample_id]

## Use the TMA Delta Distribution
delta.thresh <- comb.deltas[,.(low=quantile(delta, .05), high=quantile(delta, .95), distr=list(ecdf(delta))),by=.(ProbeName, type)]

#For each specimen, need to only include data from one batch. In other words, specimens tested in replicates will give errors, when tested across multiple dates.
#This prints out all the possible combinations
#will join with datatables that allow to see which of these are from the same patients
my.scores %>% dplyr::select(sample_id, num_batch) %>% distinct() %>% group_by(sample_id) %>% summarize(batches_per_specimen = n()) %>%
  right_join(., my.scores, by="sample_id") %>% dplyr::select(sample_id, num_batch, batches_per_specimen) %>% as.data.frame() %>% distinct() %>%
  arrange(batches_per_specimen, sample_id)

#For the filter option, one can select the specimen/batch combos
t9.dt <- my.scores %>%
  filter((sample_id == my_samp & num_batch == runid) | sample_id == my_samp_compare & num_batch == runid_compare) %>%
  rename("avg_abund" = "norm") %>% mutate(use_roi = "avg") %>%
  dplyr::select("Segment (Name/ Label)", "segment_label", "sample_id", "use_roi", "ProbeName", "avg_abund", "num_batch") %>%
  as.data.table()

#Changed color from segment # to segment_label
t9.cast <- dcast(`segment_label`+ProbeName~sample_id, value.var="avg_abund",data=t9.dt[use_roi == "avg"])
colnames(t9.cast)<- c('segment_label','ProbeName','Specimen_ID','Specimen_ID_compare')
t9.cast[,delta:=t9.cast$Specimen_ID-t9.cast$Specimen_ID_compare]
t9.cast <- merge(t9.cast, delta.thresh, by="ProbeName",allow.cartesian=T)
t9.cast[,`:=`(status="neither")]
t9.cast[delta < low, status:="low"]
t9.cast[delta > high, status:="high"]
t9.cast[,perc:=mapply(function(d, x) d(x), distr, delta)]
t9.cd <- comb.deltas[ProbeName %in% t9.cast[,ProbeName]]
t9.cast[,fix_status:=ifelse(status=="neither", "neither", "high/low")]

#Stop script if specimen(s) do not have both stromal tumor segments
stopifnot(!any(is.na(t9.cast)))

# WRITE TO PDF

# Produce Cover Sheet
tt1 <- ttheme_minimal(core=list(fg_params=list(fontface=3, fontsize=23)))
tt_cover <- ttheme_minimal(core=list(bg_params = list(fill = blues9[1:4], col=NA),
                                     fg_params=list(fontface=3, fontsize=23)),
                           colhead=list(fg_params=list(col="darkblue", fontface=4L, fontsize=30)))

summ_df <- as.data.frame(c(my_samp, my_samp_compare, runid, runid_compare, as.character(Sys.Date())), header=FALSE)
rownames(summ_df) <- c('SAMPLE ID: ', 'COMPARATOR SAMPLE ID: ', 'RUN ID: ', 'COMPARATOR RUN ID: ','RUN DATE: ')
colnames(summ_df) <- c('Nanostring_DSP')

gc1 <- tableGrob(summ_df, theme=tt_cover)
grid.draw(gc1)

#Data table with results for each Probe
# This stuff is for the table output, but we'll do it up here so the delta plots can also access the score_tum tables.
tt <- ttheme_default(base_size = 12)
score_out <- t9.cast %>% dplyr::select("ProbeName", "segment_label")
score_out <- cbind(score_out, t9.cast$Specimen_ID, t9.cast$Specimen_ID_compare)
score_out <- cbind(score_out, t9.cast$delta, t9.cast$low, t9.cast$status)
colnames(score_out) <- c("ProbeName", "segment_label", my_samp, my_samp_compare, "delta", "low", "status")

score_out<- score_out %>% mutate(across(3:6, round, 2)) %>%
  rename("Region"="segment_label", "Protein"="ProbeName", "Delta"="delta", "Threshold"="low", "Interpretation"="status") %>%
  mutate(Threshold = abs(Threshold)) %>%
  mutate(Interpretation = ifelse(Interpretation == "neither", "indeterminate", Interpretation))
score_tum <- score_out[`Region` == 'tumor']
score_str <- score_out[`Region` == 'stroma']

#Changed color from segment # to segment_label
t9.cd.pg1 <- t9.cd[`ProbeName` %in% score_tum[1:17]$Protein]
t9.cd.pg2 <- t9.cd[`ProbeName` %in% score_tum[18:34]$Protein]
t9.cd.pg3 <- t9.cd[`ProbeName` %in% score_tum[35:51]$Protein]
t9.cd.pg4 <- t9.cd[`ProbeName` %in% score_tum[52:68]$Protein]
t9.cast.pg1 <- t9.cast[`ProbeName` %in% score_tum[1:17]$Protein]
t9.cast.pg2 <- t9.cast[`ProbeName` %in% score_tum[18:34]$Protein]
t9.cast.pg3 <- t9.cast[`ProbeName` %in% score_tum[35:51]$Protein]
t9.cast.pg4 <- t9.cast[`ProbeName` %in% score_tum[52:68]$Protein]

# Plot the first page of deltas
t9.plot <- ggplot(data=t9.cd.pg1, mapping=aes(x=delta)) +
  geom_density() +
  geom_vline(data=t9.cast.pg1, mapping=aes(xintercept=delta, color=`segment_label`, linetype=fix_status)) +
  facet_grid(ProbeName~type) +
  scale_linetype_manual(values=c(neither="dashed", `high/low`="solid")) +
  geom_text_repel(data=t9.cast.pg1, mapping=aes(y=1, x=delta, label=paste0(round(delta, 2), " (", round(perc, 2), ")"))) +
  theme_bw() + ylab("") + xlab("Delta (Bx2 - Bx1)") +
  theme(strip.text.y=element_text(angle=0), axis.text.y=element_blank(), axis.ticks.y=element_blank())
grid.draw(t9.plot)

# Plot the second page of deltas
t9.plot <- ggplot(data=t9.cd.pg2, mapping=aes(x=delta)) +
  geom_density() +
  geom_vline(data=t9.cast.pg2, mapping=aes(xintercept=delta, color=`segment_label`, linetype=fix_status)) +
  facet_grid(ProbeName~type) +
  scale_linetype_manual(values=c(neither="dashed", `high/low`="solid")) +
  geom_text_repel(data=t9.cast.pg2, mapping=aes(y=1, x=delta, label=paste0(round(delta, 2), " (", round(perc, 2), ")"))) +
  theme_bw() + ylab("") + xlab("Delta (Bx2 - Bx1)") +
  theme(strip.text.y=element_text(angle=0), axis.text.y=element_blank(), axis.ticks.y=element_blank())
grid.draw(t9.plot)

t9.plot <- ggplot(data=t9.cd.pg3, mapping=aes(x=delta)) +
  geom_density() +
  geom_vline(data=t9.cast.pg3, mapping=aes(xintercept=delta, color=`segment_label`, linetype=fix_status)) +
  facet_grid(ProbeName~type) +
  scale_linetype_manual(values=c(neither="dashed", `high/low`="solid")) +
  geom_text_repel(data=t9.cast.pg3, mapping=aes(y=1, x=delta, label=paste0(round(delta, 2), " (", round(perc, 2), ")"))) +
  theme_bw() + ylab("") + xlab("Delta (Bx2 - Bx1)") +
  theme(strip.text.y=element_text(angle=0), axis.text.y=element_blank(), axis.ticks.y=element_blank())
grid.draw(t9.plot)

t9.plot <- ggplot(data=t9.cd.pg4, mapping=aes(x=delta)) +
  geom_density() +
  geom_vline(data=t9.cast.pg4, mapping=aes(xintercept=delta, color=`segment_label`, linetype=fix_status)) +
  facet_grid(ProbeName~type) +
  scale_linetype_manual(values=c(neither="dashed", `high/low`="solid")) +
  geom_text_repel(data=t9.cast.pg4, mapping=aes(y=1, x=delta, label=paste0(round(delta, 2), " (", round(perc, 2), ")"))) +
  theme_bw() + ylab("") + xlab("Delta (Bx2 - Bx1)") +
  theme(strip.text.y=element_text(angle=0), axis.text.y=element_blank(), axis.ticks.y=element_blank())
grid.draw(t9.plot)

# Plots the delta tables
g <- tableGrob(score_tum[1:34,1:7], rows = NULL, theme = tt)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 1, l = 1, r = ncol(g))
grid.newpage()
grid.draw(g)

g <- tableGrob(score_tum[35:68,1:7], rows = NULL, theme = tt)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 1, l = 1, r = ncol(g))
grid.newpage()
grid.draw(g)

g <- tableGrob(score_str[1:34,1:7], rows = NULL, theme = tt)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 1, l = 1, r = ncol(g))
grid.newpage()
grid.draw(g)

g <- tableGrob(score_str[35:68,1:7], rows = NULL, theme = tt)
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                     t = 1, l = 1, r = ncol(g))
grid.newpage()
grid.draw(g)

dev.off()
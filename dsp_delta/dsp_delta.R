# VERSION: 1.0.1
# Version history
# 1.0.1 - Named all arguments, reformatting, allow for reference to contain only segment 1

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
selected_pos <- args[6]
pos_cntrls <- args[7]

# dsp_runner intermediates
exp.out <- args[8]
seg.proc.out <- args[9]

# Output files
delta.ridge.out <- args[10]
plot.out <- args[11]

paths <- data.table(read.xlsx(ab_info, sheet="parsed"))
load(exp.out)
load(seg.proc.out)

## Pairwise Delta Distribution from RUV Normalized TMA
tma.meta[,tmpval:=1]
cl.repmat <- reshape2::acast(barcode~name, value.var="tmpval", data=tma.meta[type=="quant"], fun.aggregate = function(x) as.integer(length(x) > 0))
#and then to line 202 to create cl.pcs
cl.pcs <- compute_factors(tma.abund, cl.repmat)

use.abs <- paths[Manuscript=='y' & ProbeName %in% c("ARID1A", "Cleaved Caspase 9") == F]

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

#Create variable for pos.cntrl.dt by going to line 415 and loading positive cell line file.
pos.cntrls <- data.table(openxlsx::read.xlsx(selected_pos))
pos.cntrl.dt <- rbindlist(lapply(split(pos.cntrls, by="Target"), function(x){
  x.cl <- strsplit(x$`Low.CV.(Positive.cell.line./.tissue)`, ', ')[[1]]
  data.table(ProbeName=x$Target, name=x.cl)
}))

pos.cntrl.dt[name == "MDA-MB-231",name:="231"]
pos.cntrl.dt[name == "MDA-MB-453",name:="453"]
pos.cntrl.dt[name == "MDA-MB-468+EGF",name:="MDA468+EGF"]
pos.cntrl.dt[name == "OVCAR8+PARPi",name:="OvCar8+PARPi"]
pos.cntrl.dt[name == "THP-1",name:="THP1"]
pos.cntrl.dt[name == "MCF-7",name:="MCF7"]
pos.cntrl.dt <- pos.cntrl.dt[ProbeName %in% intersect(rownames(all.norm), use.abs[pathway %in% c("Negative Controls", "Expression Controls")==F,ProbeName])]

melt.dsp <- merge(melt.dsp, pos.cntrl.dt, by=c("ProbeName", "name"))


by.ab <- split(melt.dsp, by="ProbeName")
val.range <- range(melt.dsp$value)

# THIS TABLE SHOWS STATS FOR NORMALIZED TMA's
cv.sum <- melt.dsp[,.(sd=sd(value), mean=mean(value), mad=mad(value),.N),by=.(ProbeName, name)]

#looks like stats are calculated for all TMA-ab combos
#But because more than one TMA per ab, MAD value for a given ab is taken from the worst of its associated tmas!!!!
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
probe.mads <- data.table(openxlsx::read.xlsx(pos_cntrls))
comb.deltas[,ab_fac:=factor(ProbeName, levels=rev(probe.mads$ProbeName), ordered=T)]

#This is the first relevant image created
delta.ridge <- ggplot(data=comb.deltas, mapping=aes(x=delta, y=ab_fac, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.975)
  ) +
  scale_fill_manual(
    name = "Quantiles", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
    labels = c("(0, 0.05]", "(0.05, 0.95]", "(0.95, 1]")
  ) +
  facet_wrap(~type) +
  theme_bw() + ylab("") + xlab("Delta")
ggsave(delta.ridge, width=7, height=7, file=paste0(delta.ridge.out))
#delta.ridge

# Form antibody scores as the quantiles relevant to reference cohort
ref_samps <- my.meta[Best_Response == "Ref",unique(sample_id)]
ref_samps <- ref_samps[!ref_samps %in% my_samp]
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
cntrls <- paths[analysis_pathway %in% c("Negative Controls", "Expression Controls")]

delta.thresh <- comb.deltas[,.(low=quantile(delta, .05), high=quantile(delta, .95), distr=list(ecdf(delta))),by=.(ProbeName, type)]
delta.thresh <- delta.thresh[ProbeName %in% cntrls$ProbeName == F]

#For each specimen, need to only include data from one batch. In other words, specimens tested in replicates will give errors, when tested across multiple dates.
#This prints out all the possible combinations
#will join with datatables that allow to see which of these are from the same patients
my.scores %>% dplyr::select(sample_id, num_batch) %>% distinct() %>% group_by(sample_id) %>% summarize(batches_per_specimen = n()) %>%
  right_join(., my.scores, by="sample_id") %>% dplyr::select(sample_id, num_batch, batches_per_specimen) %>% as.data.frame() %>% distinct() %>%
  arrange(batches_per_specimen, sample_id)

#For the filter option, one can select the specimen/batch combos
t9.dt <- my.scores %>%
  #filter((sample_id == "DN21-00153B" & num_batch !="10052021") | sample_id == "DN21-00163B") %>%
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
t9.delta <- t9.cast[status != "neither"]
t9.cd <- comb.deltas[ProbeName %in% t9.delta[,ProbeName]]
t9.cast[,fix_status:=ifelse(status=="neither", "neither", "high/low")]

#Stop script if specimen(s) do not have both stromal tumor segments
stopifnot(!any(is.na(t9.cast)))

#Changed color from segment # to segment_label
t9.plot <- ggplot(data=t9.cd, mapping=aes(x=delta)) +
  geom_density() +
  geom_vline(data=t9.cast[ProbeName %in% t9.delta[,ProbeName]], mapping=aes(xintercept=delta, color=`segment_label`, linetype=fix_status)) +
  facet_grid(ProbeName~type) +
  scale_linetype_manual(values=c(neither="dashed", `high/low`="solid")) +
  geom_text_repel(data=t9.cast[ProbeName %in% t9.delta[,ProbeName]], mapping=aes(y=1, x=delta, label=paste0(round(delta, 2), " (", round(perc, 2), ")"))) +
  theme_bw() + ylab("") + xlab("Delta (Bx2 - Bx1)") +
  theme(strip.text.y=element_text(angle=0), axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave(t9.plot, file=paste0(plot.out), width=14, height=8)
t9.plot
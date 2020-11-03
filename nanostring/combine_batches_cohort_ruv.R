#!/usr/bin/env Rscript

## Takes 2 or more samples run RUV. Performs QC and plots boxplot and binmaps
## Usage ./combined_batches_ruv.R --help
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ruv))
suppressPackageStartupMessages(library(gridExtra))

version="4.0"
mat_version = "20200512"
ihc_version = "20200507"

parser <- ArgumentParser()
parser$add_argument("--data_dir", type="character",
default="/Volumes/OHSU/CLINICAL/Nanostring/Assay_Whole_Slide/output", 
                    dest="data_dir", help="directory name")
parser$add_argument("--comb_sheet", type="character", dest="comb_sheet", help="ruv samplesheet")
parser$add_argument("--validation_file",
type="character",default=paste0("/Volumes/OHSU/CLINICAL/Nanostring/Assay_Whole_Slide/REFERENCE_FILES/validation_samples_rawdata_", mat_version, ".txt"),
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--md_file", type="character", default= "/Users/patterja/Box/NANOSTRING/nanostring_metadata.xlsx",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--ihc_file", type="character",
default=paste0("/Volumes/OHSU/CLINICAL/Nanostring/Assay_Whole_Slide/REFERENCE_FILES/ihc_status_", ihc_version, ".txt"),
                    dest="ihc_file", help="ihc file")
parser$add_argument("--ab_ref_file", type="character", default="/Volumes/OHSU/CLINICAL/Nanostring/Assay_Whole_Slide/REFERENCE_FILES/ANTIBODY_REFERENCE_v1.0.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--include_ctrls", action="store_true", default=FALSE,
                    dest="include_ctrls", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
data_dir = args$data_dir
comb_sheet = args$comb_sheet
validation_file = args$validation_file
md_file = args$mbc_md_file
ab_ref_file = args$ab_ref_file
include_ctrls = args$include_ctrls
ihc_file = args$ihc_file

#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468.control","MDA468+EGF")
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2")
omitregex = paste0(paste0("^", omit), collapse = "|")
ab.ctrl = "IgG|POS|NEG|^S6|^Histone"



#FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
pcaplot <- function (mat, title = "PCA Plot", col=rownames(mat)) {
  #' pca
  #'
  #' @param mat (matrix/dataframe):  mat
  #' @param title (character) :
  #' @param subtype (dataframe metadata) : rownames correspnond to
  #' @param labe (character) : rownames correspnond to
  #' @export
  #' @return ggplot
  col = c(col)
  var = mat[apply(mat, 1, var, na.rm = TRUE) != 0, ]
  cc.var = var[complete.cases(var), ]
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  PC1_and_PC2 = data.frame(PC1 = pca_prcomp$x[, 1], PC2 = pca_prcomp$x[,2], type = rownames(pca_prcomp$x))
  perc = (pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) * 100
  labs <- sapply(seq_along(perc), function(i) {
    paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  
  PCsmd = cbind(PC1_and_PC2, col=col)
  levs = levels(factor(col))
  cols =c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#8DD3C7","#FFFFB3",
          "#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F",
          "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#8DD3C7","#FFFFB3",
          "#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
  p = ggplot(PCsmd,aes_string("PC1", "PC2", col="col")) +
    geom_point(size = 1.5) +
    geom_text(aes(label = PCsmd$type), vjust = -1, size=2) +
    labs(title = title,x = labs[1], y = labs[2]) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "gray90"),
          panel.border = element_rect(colour = "gray90", fill=NA),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 4), legend.position = "right") +
    scale_colour_manual(values =cols[1:length(levs)]) +
    xlim(-20,20) +ylim(-9,9)
  return(p)
}

rel_log_exp <- function (Y)
  #' relative log expression calculation
  #'
  #' @param Y data matrix. Rows are observations and columns are features (e.g. genes).
  #' @return numeric
  #' @export
{
  features = colnames(Y)
  print(paste0("Using columns as features.\nFeatures are: ", paste0(features[1:10], collapse=","), "..."))
  
  Yrle = apply(Y, 2, function(x) x-median(x))
  
  return(Yrle)
}
## OUTPUT PREP ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#outdirectory. Not sure how this will integrate with galaxy
outdir = dirname(normalizePath(comb_sheet))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~ samplesheet
print(paste0("output files will be in current dir: ", getwd()))

samps2batchcorr = read.table(file= comb_sheet, sep="\t", stringsAsFactors = F, header=T)
sampname = gsub("_samplesheet.txt", "", basename(comb_sheet))

#~ metadata
md = read.xlsx(xlsxFile=md_file, sheet="nanostring_metadata", check.names=T)
md$sampcolumn = make.names(paste0(md$Batch, "__", md$Sample.Name))

# antibody metadata 
ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)

#~ validation
validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = T)

#ihc status 
ihc = read.csv(ihc_file, sep= "\t", check.names = T)
ihc$sampcolumn = make.names(paste0(ihc$Batch, "__", ihc$Sample.Name))

#remove ihc duplicates and make cohorts samples ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
root_name = sapply(strsplit(as.character(ihc$Sample.Name), "\\s+"), `[`, 1)
idx_dup = duplicated(root_name)
ihc = ihc[!idx_dup,]

cohs = list(BC_preTx = ihc$sampcolumn[ihc$cohort=="breast" & ihc$TNBC=="FALSE"], 
            TNBC_preTx = ihc$sampcolumn[ihc$cohort=="breast" & ihc$TNBC=="TRUE"],
            TNBC_onTx = ihc$sampcolumn[ihc$cohort=="breast_onTx" & ihc$TNBC=="TRUE"])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ PARSE ALL SAMPLES FROM ALL BATCHES IN SAMPSHEET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
raw_dats = c()
raw_samps = c()
for (f in 1:length(unique(samps2batchcorr$Batch))){
  
  #get raw data
  batch_name = samps2batchcorr$Batch[f]
  file_name = as.character(paste0(data_dir,"/",batch_name ,"/rawdata.txt"))
  if (file.exists(file_name)){
    batch = read.table(file = file_name, sep="\t", row.names=2, stringsAsFactors=F, header=T, check.names = T)
  } else {
    
    stop(sprintf("Batch QC failed. No rawdata.txt file exists! Batch %s", batch_name))
  }
  batch[,c("CodeClass", "Accession")] <- NULL
  
  ## TEST ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #TEST to see if sample names in batch directories
  if(any(colnames(batch) %in% make.names(samps2batchcorr$Sample.Name))){
    print(paste0("Sample:", samps2batchcorr$Sample.Name[f], " in ", file_name))
    } else {
    print(paste0("Sample:", samps2batchcorr$Sample.Name[f], " NOT in ", file_name, "\n Check sample name."))
    stop()
    }
  #TEST to make sure the controls are in each batch
  for (ctrl in make.names(controls)){
    if (!ctrl %in% colnames(batch)){
      print(paste("WARNING: There is a control missing: ", ctrl))
    } else{
      print(paste0(ctrl, " controls accounted for."))
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #get only samples you want to batch and the controls
  idx_controls = which(make.names(colnames(batch)) %in% make.names(controls))
  idx_sample = which(make.names(colnames(batch)) %in% make.names(samps2batchcorr$Sample.Name))
  sel_batch = batch[,c(idx_controls, idx_sample)]
  sel_md = data.frame("batch"=batch_name,
                      "sample" = colnames(sel_batch),
                      "fullname" = make.names(paste0(batch_name,"__", colnames(sel_batch))), stringsAsFactors = F)
  
  colnames(sel_batch) = make.names(paste0(batch_name,"__", colnames(sel_batch)))
  
  # combine with all batches if not first one from samplesheet
  if (f==1){
    raw_dats = data.matrix(sel_batch)
    raw_samps = sel_md
  } else{
    raw_dats = cbind(raw_dats, sel_batch[match(rownames(raw_dats), rownames(sel_batch)),])
    raw_samps = rbind(raw_samps, sel_md)
  }
}
colnames(raw_dats) = paste0("combining_", colnames(raw_dats))
raw_samps$fullname = paste0("combining_", raw_samps$fullname)


# COMBINING DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## combine for replicate structure: rows combined with validation, ctrls and samples
combined_raw = cbind(validation, raw_dats[match(rownames(validation), rownames(raw_dats)),])
lcombined = log(combined_raw+1)

# COMBINED DATA METADATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create a new metadata for replicate structure
combined_md = data.frame(batch = sapply(strsplit(as.character(colnames(lcombined)), "__"), `[`, 1), 
                         samp =sapply(strsplit(as.character(colnames(lcombined)), "__"), `[`, 2),
                         sampcolumn = c(colnames(lcombined)),stringsAsFactors = F)
md = md[md$sampcolumn %in% combined_md$sampcolumn,]
combined_md$deid = ifelse(combined_md$samp %in% make.names(controls), yes=combined_md$samp, no="sample")



#~ RUV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
combined_md$reps = ifelse(combined_md$samp %in% make.names(controls), yes=combined_md$samp, no=combined_md$sampcolumn)
repmat = replicate.matrix(combined_md$reps)
idx_ctl = grep(ab.ctrl, rownames(lcombined))

RUVcorrected = RUVIII(Y=t(lcombined), ctl=idx_ctl, M=repmat, k=2, include.intercept = FALSE, inputcheck = F)
bcdat = t(RUVcorrected[,!grepl("POS|NEG", rownames(lcombined))])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## because i got sick of calling it all the time
celllines = combined_md$sampcolumn[combined_md$samp %in% make.names(controls)]
valid_celllines = setdiff(celllines, raw_samps$fullname)
samps = raw_samps$fullname[!raw_samps$sample %in% make.names(controls)]
samp_celllines = setdiff(raw_samps$fullname, samps)


#~ QC PLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("ruv_figures", showWarnings = F)

#~ RLE box plot
lctls=lcombined[!grepl("POS|NEG", rownames(lcombined)),celllines]

limits_rle_raw = max(as.matrix(lctls))-median(as.matrix(lctls))
rle_rawbox = ruv_rle(Y = t(lctls), rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% celllines,]), 
                  ylim=c(-limits_rle_raw,limits_rle_raw)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  geom_text(aes(x = rle.x.factor, y=-limits_rle_raw-0.5, label=colnames(lctls)), angle=90, hjust=0, size=2)+
  theme(legend.position = "right", legend.text = element_text(size=6)) + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression \nNormalized Data (IGG corrected) \n", sampname, "_", "PRE-RUV"))

ctl_bc = bcdat[,celllines]
limits_rle_bc = max(as.matrix(ctl_bc))-median(as.matrix(ctl_bc))

rle_ruvbox = ruv_rle(Y = t(ctl_bc), 
                   rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% celllines,]), 
                   ylim=c(-limits_rle_bc,limits_rle_bc)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  geom_text(aes(x = rle.x.factor, y=-limits_rle_bc-0.5, label="RUVcorrected"), angle=90, hjust=0, size=2) +
  theme(legend.position = "right") + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression RUV processed\n", sampname, "_", "RUVcorrected"))



# NORMAL HEATMAP
pdf(paste0("ruv_figures/",sampname,"_control_heatmap.pdf"), width = 7, height = 8)
## rle raw
rownames(combined_md) = combined_md$sampcolumn
pheat_raw = pheatmap(mat = t(lctls),
                     color             = colorRampPalette(c("green", "black", "red"))(50),
                     border_color      = NA,
                     show_colnames     = TRUE,
                     show_rownames     = TRUE,
                     annotation_row    = combined_md[colnames(lctls),c("reps"), drop=F], 
                     drop_levels       = TRUE,
                     fontsize          = 8,
                     fontsize_row      = 5,
                     fontsize_col      = 5,
                     main              = "log raw", 
                     cluster_rows      = T,
                     cluster_cols      = T)

## ruv raw
pheat_ruv = pheatmap(mat = t(ctl_bc),
                     color             = colorRampPalette(c("green", "black", "red"))(50),
                     border_color      = NA,
                     show_colnames     = TRUE,
                     show_rownames     = TRUE,
                     annotation_row    = combined_md[colnames(lctls),c("reps"), drop=F], 
                     drop_levels       = TRUE,
                     fontsize          = 8,
                     fontsize_row      = 5,
                     fontsize_col      = 5,
                     main              = "RUV corrected controls", 
                     cluster_rows      = T,
                     cluster_cols      = T)
dev.off()
# RLE HEATMAP
## rle raw
pdf(paste0("ruv_figures/",sampname,"_RLE_heatmap.pdf"), width = 7, height = 8)
rleraw = rel_log_exp(t(lctls))
pheat_rleraw = pheatmap(
  mat = rleraw,
  color             = colorRampPalette(c("green", "black", "red"))(50),
 border_color      = NA,
 show_colnames     = TRUE,
 show_rownames     = TRUE,
 annotation_row    = combined_md[colnames(lctls),c("reps"), drop=F], 
 drop_levels       = TRUE,
 fontsize          = 8,
 fontsize_row      = 5,
 fontsize_col      = 5,
 main              = "rle raw", 
 cluster_rows      = T,
 cluster_cols      = T)

## ruv rle
rleruv = rel_log_exp(t(ctl_bc))

pheat_rleruv = pheatmap(mat = (rleruv),
                     color             = colorRampPalette(c("green", "black", "red"))(50),
                     border_color      = NA,
                     show_colnames     = TRUE,
                     show_rownames     = TRUE,
                     annotation_row    = combined_md[colnames(lctls),c("reps"), drop=F], 
                     drop_levels       = TRUE,
                     fontsize          = 8,
                     fontsize_row      = 5,
                     fontsize_col      = 5,
                     main              = "RLE of RUV corrected controls", 
                     cluster_rows      = T,
                     cluster_cols      = T,
                     filename = )
dev.off()
#~ PCA plots
pcacol= combined_md$deid[match(colnames(lctls), combined_md$sampcolumn)]
pca_raw = pcaplot(mat = lctls, title = paste0("log2 (Normalized Signal + 1)\n", sampname, "_", "PRE-RUV"), col=pcacol)

pcacol= combined_md$deid[match(colnames(ctl_bc), combined_md$sampcolumn)]
pca_ruv = pcaplot(ctl_bc, title = paste0("RUV Processed\n", sampname, "_", "RUVcorrected"), pcacol)

#~ TRA 
#t(RUVcorrected) is post norm, comb is pre norm, using no antibody filtered data
eruv = t(exp(ctl_bc))
rawctls = t(combined_raw)[rownames(eruv),colnames(eruv)]
tra = matrix(NA, nrow = 0, ncol = 4)

for (ctrl in make.names(controls)){
  #get matrix of celllines for controls not including batches of interest
  valid.ctrl_names = combined_md$sampcolumn[combined_md$samp == ctrl]

  selv.ctrl_raw = rawctls[valid.ctrl_names,, drop=F]
  selv.ctrl_ruv = eruv[valid.ctrl_names,, drop=F]
  tra_ctrlraw =c()
  tra_ctrlruv =c()
  for (r in 1:nrow(selv.ctrl_ruv)){
    tra_ctrlraw = rbind(tra_ctrlraw, log(sweep(as.matrix(selv.ctrl_raw[-r,,drop=F]), 2, as.numeric(selv.ctrl_raw[r,,drop=F]), `/`)))
    tra_ctrlruv = rbind(tra_ctrlruv, log(sweep(as.matrix(selv.ctrl_ruv[-r,,drop=F]), 2, as.numeric(selv.ctrl_ruv[r,,drop=F]), `/`)))
  }
  tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrlraw)), "sample"="raw"))
  tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrlruv)), "sample"="ruv"))
}

tra_btwn_batch = matrix(NA, nrow = 0, ncol = 4)

allselraw = t(combined_raw[!grepl("POS|NEG", rownames(lcombined)),samp_celllines])
allselruv = t(exp(bcdat[,samp_celllines]))

for (ctrl in make.names(controls)){
  #get matrix of celllines for controls not including batches of interest
  valid.ctrl_names = samp_celllines[grep(ctrl, samp_celllines)]
  if (length(valid.ctrl_names)==1){
    selv.ctrl_raw = rbind(allselraw[valid.ctrl_names,, drop=F],allselraw[valid.ctrl_names,, drop=F])
    selv.ctrl_ruv = rbind(allselruv[valid.ctrl_names,, drop=F],allselruv[valid.ctrl_names,, drop=F])
  } else{
    selv.ctrl_raw = allselraw[valid.ctrl_names,, drop=F]
    selv.ctrl_ruv = allselruv[valid.ctrl_names,, drop=F]
  }
  tra_ctrlraw =c()
  tra_ctrlruv =c()
  for (r in 1:nrow(selv.ctrl_ruv)){
    tra_ctrlraw = rbind(tra_ctrlraw, log(sweep(as.matrix(selv.ctrl_raw[-r,,drop=F]), 2, as.numeric(selv.ctrl_raw[r,,drop=F]), `/`)))
    tra_ctrlruv = rbind(tra_ctrlruv, log(sweep(as.matrix(selv.ctrl_ruv[-r,,drop=F]), 2, as.numeric(selv.ctrl_ruv[r,,drop=F]), `/`)))
  }
  tra_btwn_batch = rbind(tra_btwn_batch, data.frame(melt(as.matrix(tra_ctrlraw)), "sample"="raw"))
  tra_btwn_batch = rbind(tra_btwn_batch, data.frame(melt(as.matrix(tra_ctrlruv)), "sample"="ruv"))
}

ptra=ggplot(tra, aes(x=sample, y=value, fill=sample)) + 
  geom_boxplot() +
  facet_wrap(~Var2, scale="free") +
  geom_hline(yintercept =0, color="red") +
  labs(title=paste0("Distribution of TRA (technical replicate agreement) RAW and RUV Corrected\nTRAs calculated for new batch control versus each of the validation replicates\n"),
       y="log(validation count/new batch count)") +
  theme(legend.position = "bottom")

ptra_btwn=ggplot(tra_btwn_batch, aes(x=sample, y=value, fill=sample)) + 
  geom_boxplot() +
  geom_hline(yintercept =0, color="red") +
  labs(title=paste0("Distribution of TRA between RAW and RUV Corrected\n between Controls of Batches being compared"),
       y="log(validation count/new batch count)") +
  theme(legend.position = "bottom")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


print("saving batch correction plots")
pdf(file=paste0("ruv_figures/",sampname,"_TRA.pdf"), width = 8, height = 8)
print(ptra)
print(ptra_btwn)
dev.off()

pdf(paste0("ruv_figures/",sampname,"_RLE_boxplots.pdf"), width = 8, height = 4)
print(rle_rawbox)
print(rle_ruvbox)
dev.off()

pdf(paste0("ruv_figures/",sampname,"_PCA.pdf"), width = 8, height = 8)
print(pca_raw)
print(pca_ruv)
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ BATCH CORRECTED DATA FILTERING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (include_ctrls){
  print("keeping all antibodies: IgG and antibodies that did not perform well or did not have dynamic range")
  fbcdat = bcdat
  other_abs=setdiff(colnames(bcdat),ab_ref$X.AbID)
  ab_order = c(ab_ref$X.AbID[order(ab_ref$Target)], other_abs)
  pathways = ab_ref$Pathway[ab_ref$X.AbID %in% ab_order]
} else {
  fbcdat = bcdat[!grepl(omitregex, rownames(bcdat)),]
  #AB_ORDER
  ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
  ab_order = ab_order[!grepl(omitregex, ab_order)]
  pathways = ab_ref$Pathway[ab_ref$X.AbID %in% ab_order]
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ SPLITTING AND MELTING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
validsamp_dat = fbcdat[,!colnames(fbcdat) %in% valid_celllines & !grepl("combining", colnames(fbcdat)),drop=F]
samps_dat = fbcdat[,!colnames(fbcdat) %in% celllines & grepl("combining", colnames(fbcdat)),drop=F]

validsamp.datm = melt(validsamp_dat,  id.vars=row.names)
samps.datm = data.frame(melt(samps_dat, id.vars=row.names))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

comb_cohorts = c()
samp_cohorts = c()
for (c in seq(1:length(cohs))){
  temp_coh = validsamp.datm[validsamp.datm$Var2 %in% make.names(cohs[[c]]),]
  print(length(table(as.character(temp_coh$Var2))))
  temp_coh$ab = names(cohs[c])
  comb_cohorts = rbind(comb_cohorts, temp_coh)
  
}
comb_cohorts$ab_Var1 = paste0(comb_cohorts$ab, "_", comb_cohorts$Var1)
samp_cohorts = samps.datm


# Get percentiles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coh_percs = list()
for (cidx in seq(1:length(unique(comb_cohorts$ab)))){
  temp = c()
  selcoh = comb_cohorts[comb_cohorts$ab==names(cohs)[cidx],]
  
  for (ab in unique(comb_cohorts$Var1)){
    selab= selcoh[selcoh$Var1==ab , "value"]
    temp = rbind(temp, ab=quantile(selab, c(0.05, 0.15, 0.25, 0.75, 0.85, 0.95)))
    rownames(temp)[which(rownames(temp)=="ab")]=ab
    
  }
  coh_percs = append(coh_percs, list(temp))
  names(coh_percs)[cidx] = names(cohs)[cidx]
}

# ANTIBODY THRESHOLDING setting with detectable flag ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## if signal for antibody below sample igg then turned to minimum of validation  cohort
## if signal above do nothing
samp_cohorts['detectable']=TRUE
for (samp in samps){
  newsamp = combined_raw[,samp, drop=F]
  #newsamp = lcombined[,samp, drop=F]
  rbigg = newsamp[which(rownames(newsamp)=="RbAb-IgG"),]
  mmigg = newsamp[which(rownames(newsamp)=="MmAb-IgG1"),]
  for (i in seq(1:length(ab_order))){
    ab = (ab_order)[i]
    idx = which(samp_cohorts$Var1==ab & samp_cohorts$Var2==samp)
    minval = min(as.numeric(comb_cohorts[comb_cohorts$Var1==ab, "value"]), na.rm=T)-1
    
    if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="rabbit"){
      val = newsamp[ab,]-rbigg
    } else if (ab_ref$Host[which(ab_ref$X.AbID==ab)]=="mouse"){
      val = newsamp[ab,]-mmigg
    } else {
      print("double check antibody name")
      stop()
    }
    if (val < 0){
      #get index of samples with ab and samp
      samp_cohorts[idx,"detectable"] = FALSE
      samp_cohorts[idx,"newvalue"] = minval
    } else{
      samp_cohorts[idx,"newvalue"] = samp_cohorts[idx,"value"]
    }
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##multiple sample thresholding  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

thresh=data.frame("5_95"=c("5%", "95%"), "15_85"=c("15%", "85%"), "25_75"=c("25%","75%"), stringsAsFactors = F)
cohorder = c()
for (cidx in seq(1:length(coh_percs))){
  coh = names(coh_percs)[cidx]
  print(coh)
  for (ab in unique(comb_cohorts$Var1)){
    #print(ab)
    cohtemp = coh_percs[[coh]][ab,]
    for (tidx in seq(1:ncol(thresh))){
      lowthresh = cohtemp[thresh[1,tidx]]
      highthresh = cohtemp[thresh[2,tidx]]
      classname = paste0(names(coh_percs)[cidx], "__", colnames(thresh)[tidx])
      for (samp in samps){
        idxab = which(samp_cohorts$Var1==ab & samp_cohorts$Var2==samp)
        samp_cohorts[idxab,classname] = ifelse(samp_cohorts[idxab,"detectable", drop=T]==FALSE, yes = "ND",
                                               no = ifelse(samp_cohorts[idxab,"value",drop=T] <= lowthresh, yes="low",
                                                           no = ifelse(samp_cohorts[idxab,"value",drop=T] >= highthresh, yes="high", no="indeterminate")))
        
        if(!classname %in% cohorder) {cohorder = c(cohorder, classname)}
      }
    }
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ BINMAPS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ab_order = intersect(ab_ref$X.AbID[order(ab_ref$Target)], samp_cohorts$Var1)

msamp = melt(samp_cohorts, id.vars = list("Var1", "Var2"), measure.vars = which(colnames(samp_cohorts) %in% cohorder))
msamp$sampnames = paste0(make.names(sapply(strsplit(as.character(msamp$Var2), "__"), `[`, 2)))
msamp$cohortnames = paste0(make.names(sapply(strsplit(as.character(msamp$variable), "__"), `[`, 1)))
msamp$names = paste0(msamp$cohortnames,"__", msamp$sampnames)


samphmlist <- list()
for (samp in samps){
  for (thresholds in colnames(thresh)){
    mmsamp = msamp[grepl(thresholds, msamp$variable) & msamp$Var2==samp,]
    mmsamp$Var1 = factor(mmsamp$Var1, levels=ab_order)
    #mmsamp$discrete_ab = as.numeric(factor(mmsamp$Var1))
    mmsamp$cohortnames = factor(mmsamp$cohortnames, levels = unique(mmsamp$cohortnames))
    hm= ggplot(mmsamp) +
      geom_tile(aes(x=cohortnames, y=Var1, fill=factor(value)),colour = "grey50") +
      scale_fill_manual(values=c("indeterminate"="white", "ND"="black","low"="yellow","high"="blue")) +
      scale_x_discrete(name ="Cohort", breaks = mmsamp$cohortnames, labels = mmsamp$names) +
      #scale_y_continuous(breaks=mmsamp$discrete_ab, label=mmsamp$Var1, sec.axis = sec_axis(~.,breaks = mmsamp$discrete_ab,labels = round(mmsamp$rat, 1), name=paste0("Ratio of Biopsies: ", rationame)))+
      labs(x="Cohort", title=paste0("Threshold on ", thresh[1,thresholds], " and ", thresh[2,thresholds],"\n",  unique(mmsamp$sampnames)), y="Antibody") +
      theme(axis.text.x = element_text(size=6, colour=c("black"), angle = 90, hjust=1, vjust=0.5),
            axis.text.y = element_text(size=6),
            axis.title.x = element_text(size=6),
            axis.title.y = element_text(size=6),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "gray60"),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size=9),
            legend.title = element_blank()) 
    
    samphmlist[[samp]][[thresholds]] = hm
  }
}
pdf(file=sprintf("%s_binmaps.pdf", sampname), width=10, height=6)
for (samp in samps){
  grid.arrange(samphmlist[[samp]][["X5_95"]], samphmlist[[samp]][["X15_85"]], samphmlist[[samp]][["X25_75"]], ncol=3, top=samp)
}

dev.off()

#~ ORDERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comb_cohorts_box = comb_cohorts
samp_cohorts_box = c()
for (coh in unique(comb_cohorts$ab)){
  samp_cohorts_box = rbind(samp_cohorts_box, data.frame(samp_cohorts, "ab"=coh))}
comb_cohorts_box$Var1 = factor(comb_cohorts$Var1, levels=ab_order)
samp_cohorts_box$Var1 = factor(samp_cohorts_box$Var1, levels=ab_order)
#~ BOXPLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~new
palette = c("#FF0000FF","#004CFFFF","#00A600FF","#984ea3","#ff7f00","#a65628")

pwboxlist = list()

for (pathway in unique(ab_ref$Pathway[order(ab_ref$Pathway.Order)])){
  if (!pathway %in% c("Non-specific")){
    ab_o = order(ab_ref$Target[ab_ref$Pathway==pathway])
    antibodies = ab_ref$X.AbID[ab_ref$Pathway==pathway][ab_o]
    comb_cohorts_box = comb_cohorts[comb_cohorts$Var1 %in% antibodies,]
    comb_cohorts_box$Var1 = factor(comb_cohorts_box$Var1, levels=antibodies)
    
    samp_cohorts_boxab = samp_cohorts_box[samp_cohorts_box$Var1 %in% antibodies,]
    samp_cohorts_boxab$Var2 = as.factor(gsub("combining_X","", samp_cohorts_boxab$Var2))
    colortable = setNames(palette[1:length(levels(samp_cohorts_boxab$Var2))],nm=levels(samp_cohorts_boxab$Var2))

    bp = ggplot(comb_cohorts_box, aes(x=ab, y=value)) + 
      geom_boxplot() +
      facet_wrap(~Var1, scale="free", nrow=1) +
      geom_point(data=samp_cohorts_boxab[samp_cohorts_boxab$detectable==TRUE,], mapping=aes(x=ab, y=newvalue, colour=Var2),shape=8, size=2, position=position_jitter(width=c(0.01)))+
      geom_text(data=samp_cohorts_boxab[samp_cohorts_boxab$detectable==FALSE,], mapping=aes(x=ab, y=newvalue, colour=Var2), label="ND", size=3, position=position_jitter(width=c(0.02))) +
      scale_color_manual(values = colortable) +
      labs(y="log batch corrected counts", title=pathway) +
      theme(panel.background = element_rect(fill = "white"),
            panel.grid.major=element_line(colour="gray90"),
            plot.title = element_text(hjust = 0.5, vjust=0, face="bold", size=10),
            legend.position="none",
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=6),
            axis.text.x = element_text(size=6, colour=c("black"),angle=15),
            axis.text.y = element_text(size=6, colour="black"),
            plot.margin=unit(c(-0.1,0.1,-0.1,0.1), "cm"))
    
    pwboxlist[[pathway]] = bp
  }
}
#grab the legend from to add it to the combined plot
legbp =  ggplot() + 
  geom_point(data=samp_cohorts_boxab[samp_cohorts_boxab$detectable==TRUE,], mapping=aes(x=ab, y=newvalue, colour=Var2),shape=8, size=2, position=position_jitter(width=c(0.01)))+
  geom_text(data=samp_cohorts_boxab[samp_cohorts_boxab$detectable==FALSE,], mapping=aes(x=ab, y=newvalue, colour=Var2), label="ND", size=3, position=position_jitter(width=c(0.02))) +
  scale_color_manual(name="BatchName__SampleName", values = colortable) +
  theme(legend.position="right")
        

pwboxlist[["legend"]]=get_legend(legbp)
lay <- rbind(c(1,NA,2,2,2,2),
             c(3,3,3,3,3,3),
             c(4,4,4,4,NA,NA),
             c(5,5,5,7,7,7),
             c(6,6,6,6,6,NA))

pdf(file=sprintf("%s_boxplots.pdf", sampname), width=10, height=10)

grid.arrange(grobs=pwboxlist,layout_matrix=lay, top=sampname)
dev.off()

#~ HEATMAP #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdf(file=paste0(sampname, "_samples_heatmap.pdf"), width = 4.5, height = 10)

pheatmap(mat = samps_dat[,samps],
         color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
         cluster_cols      = T,
         cluster_rows      = T,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         main              = sampname,
         labels_col = gsub("Batch Corrected combining_", "", colnames(samps_dat[,samps])))


dev.off()


#~ SCATTER PAIRS 
#colnames(plt_samps) = gsub("^.*1020_", "", colnames(plt_samps))

pdf(file=paste0(sampname, "_sample_correlation.pdf"), width = 8, height = 6)

pairs(samps_dat, lower.panel = function(x,y,...){points(x,y,...);abline(a = 0,b = 1,...)}, upper.panel = function(x, y) {usr <- par("usr"); on.exit(par(usr)); par(usr = c(0, 1, 0, 1)); r <- round(cor(x, y), digits=2); txt <- paste0("R = ", r); text(0.5, 0.5, txt)},main=paste0("RUV corrected Correlation of ", sampname))

print("done with comparison plots saving tables")

#~ TABLES #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(x = t(bcdat), file=file.path("ruv_figures", paste0(sampname, "_RUVcorrected.txt")),sep="\t", quote = F, row.names = T, col.names = NA)
write.table(x = t(combined_raw), file=file.path("ruv_figures", paste0(sampname, "_raw.txt")), sep="\t", quote = F, row.names = T, col.names = NA)


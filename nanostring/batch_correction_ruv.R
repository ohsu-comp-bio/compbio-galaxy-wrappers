#!/usr/bin/env Rscript

## RUV batch correction
## 
## Usage: ./ruv_batchcorrection.R --help 
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ruv))
suppressPackageStartupMessages(library(openxlsx))
#suppressPackageStartupMessages(library(nanoprot))
suppressPackageStartupMessages(library(gridExtra))

## ARGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
version="4.0"
mat_version = "20200512"
ihc_version = "20200507"
parser <- ArgumentParser()

parser$add_argument("--input", help="rawdata.txt", dest="input_file")
parser$add_argument("--validation_file", type="character", 
                    default=paste0("/Volumes/OHSU/CLINICAL/Nanostring/Assay_Whole_Slide/REFERENCE_FILES/validation_samples_rawdata_", mat_version, ".txt"),
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--md_file", type="character", default= "/Users/patterja/Box/NANOSTRING/nanostring_metadata.xlsx",
                    dest="md_file", help="metadata file")
parser$add_argument("--ihc_file", type="character",
default=paste0("/Volumes/OHSU/CLINICAL/Nanostring/Assay_Whole_Slide/REFERENCE_FILES/ihc_status_", ihc_version, ".txt"),
                    dest="ihc_file", help="ihc file")
parser$add_argument("--ab_ref_file", type="character", default=
"/Volumes/OHSU/CLINICAL/Nanostring/Assay_Whole_Slide/REFERENCE_FILES/ANTIBODY_REFERENCE_v1.0.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--include_ctrls", action="store_true", default=FALSE,
                    dest="include_ctrls", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
input_file = args$input_file
validation_file = args$validation_file
md_file = args$md_file
ab_ref_file = args$ab_ref_file
include_ctrls = args$include_ctrls
ihc_file = args$ihc_file

## FUNCTIONS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

dir.create("ruv_figures", showWarnings = F)

## INPUT DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# metadata
md = read.xlsx(xlsxFile=md_file, sheet="nanostring_metadata", check.names=T)
md$sampcolumn = make.names(paste0(md$Batch, "__", md$Sample.Name))

# validation file
validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = T)

#ihc status
ihc = read.csv(ihc_file, sep= "\t", check.names = T)
ihc$sampcolumn = make.names(paste0(ihc$Batch, "__", ihc$Sample.Name))

#remove ihc duplicates
root_name = sapply(strsplit(as.character(ihc$Sample.Name), "\\s+"), `[`, 1)
idx_dup = duplicated(root_name)
ihc = ihc[!idx_dup,]

cohs = list(BC_preTx = ihc$sampcolumn[ihc$cohort=="breast" & ihc$TNBC=="FALSE"], 
            TNBC_preTx = ihc$sampcolumn[ihc$cohort=="breast" & ihc$TNBC=="TRUE"],
            TNBC_onTx = ihc$sampcolumn[ihc$cohort=="breast_onTx" & ihc$TNBC=="TRUE"])

# antibody metadata
ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)
pathways = data.frame(ab_ref$X.AbID, Pathway = sapply(strsplit(as.character(ab_ref$Pathway), ","), `[`, 1))
pathways = rbind(pathways, data.frame(ab_ref$X.AbID, Pathway=sapply(strsplit(as.character(ab_ref$Pathway), ","), `[`, 2)))
pathways = rbind(pathways, data.frame(ab_ref$X.AbID, Pathway=sapply(strsplit(as.character(ab_ref$Pathway), ","), `[`, 3)))

pathways = pathways[complete.cases(pathways),]

# NEW BATCH
new_batch = read.table(file = input_file, sep="\t", row.names=2, stringsAsFactors=F, header=T, check.names = T)
new_batch[,c("CodeClass", "Accession")] <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GLOBAL things to include and exclude
controls=make.names(c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF"))
ab.ctrl = make.names(rownames(new_batch)[grepl("IgG|POS|NEG|^S6|^Histone", rownames(new_batch))])
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2", "NEG", "POS")
omitregex = paste0(paste0("^", omit), collapse = "|")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOOP THRU EACH SAMP #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (samp in setdiff(colnames(new_batch), controls)) {
  print(samp)
  
  # COMBINING: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # rows combined with validation  cohort, columns=ctrls and 1 samp only
  idx_controls = which(colnames(new_batch) %in% controls)
  newsampctl = new_batch[,c(idx_controls,which(colnames(new_batch)==samp))]
  colnames(newsampctl) = paste0("newbatch", "__", colnames(newsampctl))
  comb = cbind(validation, newsampctl[match(rownames(validation), rownames(newsampctl)),])

  # METADATA and RUV: adjust metadata to match ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  comb_md = data.frame(batch = make.names(sapply(strsplit(as.character(colnames(comb)), "__"), `[`, 1)), 
                       samp =make.names(sapply(strsplit(as.character(colnames(comb)), "__"), `[`, 2)),
                       sampcolumn = c(colnames(comb)),stringsAsFactors = F, check.names = T)
  md = md[md$sampcolumn %in% comb_md$sampcolumn,]
  valid_controls=make.names(c(md$sampcolumn[md$Study=="control"],paste0("newbatch__", controls)))
  
  # REPLICATE STRUCTURE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # repmat based only on controls in validation cohort for RUV input
  comb_md$reps = ifelse(comb_md$sampcolumn %in% valid_controls, yes=comb_md$samp, no=comb_md$sampcolumn)

  #LOG 
  lcomb=log(comb+1)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ BATCH CORRECTION: RUV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  repmat = replicate.matrix(comb_md$reps)
  negctl = which(make.names(rownames(lcomb)) %in% ab.ctrl)
  # r by c matrix, r=observations, c=features
  RUVcorrected = RUVIII(Y=t(lcomb), ctl=negctl, M=repmat, k=2, include.intercept = FALSE, inputcheck = FALSE)
  ruvmat = RUVcorrected[,!grepl("NEG|POS", colnames(RUVcorrected))]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ PLOTTING QC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LOG RAW   ## filter out negatives and ab_order
  lctls = lcomb[!grepl("NEG|POS", rownames(lcomb)),] 
  lctls = lctls[!make.names(rownames(lctls)) %in% ab.ctrl, comb_md$sampcolumn[comb_md$samp %in% controls]]
  limits_rle_raw = max(as.matrix(lctls))-median(as.matrix(lctls))
  # RUV 
  ctl_ruv = t(ruvmat)[, comb_md$sampcolumn[comb_md$samp %in% controls]]
  limits_rle_ruv = max(as.matrix(ctl_ruv))-median(as.matrix(ctl_ruv))
  
  # RLE BOXPLOTS
  ## rle raw
  rle_raw = ruv_rle(Y = t(lctls), 
                     rowinfo = as.matrix(comb_md[comb_md$samp %in% controls,]), 
                     ylim=c(-limits_rle_raw,limits_rle_raw)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle_raw-0.5, label=samp), angle=90, hjust=0, size=2)+
    theme(legend.position = "right", legend.text = element_text(size=6)) + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression Raw Data",  samp))
  
  ## rle batch corrected
  rle_ruv = ruv_rle(Y = t(ctl_ruv[!make.names(rownames(lctls)) %in% ab.ctrl,]), 
                    rowinfo = as.matrix(comb_md[comb_md$samp %in% controls,]), 
                    ylim=c(-limits_rle_ruv,limits_rle_ruv)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle_ruv-0.5, label=samp), angle=90, hjust=0, size=2) +
    theme(legend.position = "right") + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression Batch Corrected\n", samp))
  
  
  # signal HEATMAP
  pdf(paste0("ruv_figures/",samp,"_control_heatmap.pdf"), width = 7, height = 8)
  rownames(comb_md) = comb_md$sampcolumn
  pheat_raw = pheatmap(mat = t(lctls),
    color             = colorRampPalette(c("green", "black", "red"))(50),
    border_color      = NA,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    annotation_row    = comb_md[colnames(lctls),c("reps"), drop=F], 
    drop_levels       = TRUE,
    fontsize          = 8,
    fontsize_row      = 5,
    fontsize_col      = 5,
    main              = "log raw", 
    cluster_rows      = T,
    cluster_cols      = T)
  
  ## NORMAL ruv
  pheat_ruv = pheatmap(mat = t(ctl_ruv),
                       color             = colorRampPalette(c("green", "black", "red"))(50),
                       border_color      = NA,
                       show_colnames     = TRUE,
                       show_rownames     = TRUE,
                       annotation_row    = comb_md[colnames(lctls),c("reps"), drop=F], 
                       drop_levels       = TRUE,
                       fontsize          = 8,
                       fontsize_row      = 5,
                       fontsize_col      = 5,
                       main              = "RUV batch corrected signal", 
                       cluster_rows      = T,
                       cluster_cols      = T)
  dev.off()
  
  # RLE HEATMAP
  pdf(paste0("ruv_figures/",samp,"_RLE_heatmap.pdf"), width = 7, height = 8)
  pheat_raw = pheatmap(mat = rel_log_exp(t(lctls)),
                       color             = colorRampPalette(c("green", "black", "red"))(50),
                       border_color      = NA,
                       show_colnames     = TRUE,
                       show_rownames     = TRUE,
                       annotation_row    = comb_md[colnames(lctls),c("reps"), drop=F], 
                       drop_levels       = TRUE,
                       fontsize          = 8,
                       fontsize_row      = 5,
                       fontsize_col      = 5,
                       main              = "RLE of log raw ", 
                       cluster_rows      = T,
                       cluster_cols      = T)
  
  ## rle ruv
  pheat_ruv = pheatmap(mat = rel_log_exp(t(ctl_ruv)),
                       color             = colorRampPalette(c("green", "black", "red"))(50),
                       border_color      = NA,
                       show_colnames     = TRUE,
                       show_rownames     = TRUE,
                       annotation_row    = comb_md[colnames(lctls),c("reps"), drop=F], 
                       drop_levels       = TRUE,
                       fontsize          = 8,
                       fontsize_row      = 5,
                       fontsize_col      = 5,
                       main              = "RLE of RUV batch corrected", 
                       cluster_rows      = T,
                       cluster_cols      = T)
  dev.off()
  
  #~ pca
  pca_raw = pcaplot(mat = lctls, 
                    title="Normalized, Pre-RUV", 
                    col = comb_md$samp[comb_md$samp %in% controls])
  pca_ruv = pcaplot(mat = ctl_ruv, 
                   title="RUV processed", 
                   col = comb_md$samp[comb_md$samp %in% controls])
 
  #~ TRA 
  #t(RUVcorrected) is post norm, comb is pre norm, using no antibody filtered data
  eruv = t(exp(ctl_ruv))
  rawctls = t(comb[colnames(eruv),rownames(eruv)])
  tra = matrix(NA, nrow = 0, ncol = 4)
  
  for (ctrl in make.names(controls)){
    #get matrix of celllines for controls not including batches of interest
    valid.ctrl_names = comb_md$sampcolumn[comb_md$samp == ctrl & !comb_md$batch=="newbatch"]
    samp.ctrl_names = comb_md$sampcolumn[comb_md$samp == ctrl & comb_md$batch=="newbatch"]
    
    selv.ctrl_raw = rawctls[valid.ctrl_names,, drop=F]
    sels.ctrl_raw = rawctls[samp.ctrl_names,, drop=F]
    selv.ctrl_ruv = (eruv[valid.ctrl_names,, drop=F])
    sels.ctrl_ruv = (eruv[samp.ctrl_names,, drop=F])
    
    tra_ctrl_raw = log(sweep(as.matrix(selv.ctrl_raw), 2, as.numeric(sels.ctrl_raw), `/`))
    tra_ctrl_ruv = log(sweep(as.matrix(selv.ctrl_ruv), 2, as.numeric(sels.ctrl_ruv), `/`))
    
    
    tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrl_raw)), "sample"="raw"))
    tra = rbind(tra, data.frame(melt(as.matrix(tra_ctrl_ruv)), "sample"="ruv"))
    }
  tra = tra[!is.infinite(tra$value) & !is.na(tra$value),]
  ptra=ggplot(tra, aes(x=sample, y=value, fill=sample)) + 
    geom_boxplot() +
    facet_wrap(~Var2, scale="free") +
    geom_hline(yintercept =0, color="red") +
    labs(title=paste0("Distribution of TRA (technical replicate agreement) RAW and RUV Corrected\nTRAs calculated for new batch control versus each of the validation replicates\n", samp, "_",ctrl),
         y="log(validation count/new batch count)") +
    theme(legend.position = "bottom")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print("saving batch correction plots")
  pdf(file=paste0("ruv_figures/",samp,"_TRA.pdf"), width = 8, height = 8)
  plot(ptra)
  dev.off()
  
  pdf(paste0("ruv_figures/",samp,"_RLE_boxplots.pdf"), width = 8, height = 4)
  print(rle_raw)
  print(rle_ruv)
  dev.off()
  
  
  pdf(paste0("ruv_figures/",samp,"_PCA.pdf"), width = 8, height = 8)
  print(pca_raw)
  print(pca_ruv)
  dev.off()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ BATCH CORRECTED DATA FILTERING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (include_ctrls){
    print("keeping all antibodies: IgG and antibodies that did not perform well or did not have dynamic range")
    nmat = t(ruvmat)
    #AB_ORDER
    other_abs=setdiff(colnames(ruvmat),ab_ref$X.AbID)
    ab_order = c(ab_ref$X.AbID[order(ab_ref$Target)], other_abs)
  } else {
    nmat = t(ruvmat)
    nmat = nmat[!grepl(omitregex, rownames(nmat)),]
    #AB_ORDER
    ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
    ab_order = ab_order[!grepl(omitregex, ab_order)]
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ SPLITTING AND MELTING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # matrix of samples only, remove controls
  val_nmat = nmat[,!colnames(nmat) %in% valid_controls & !grepl("newbatch", colnames(nmat)),drop=F]
  new_nmat = nmat[,!colnames(nmat) %in% valid_controls & grepl("newbatch", colnames(nmat)),drop=F]
  #~ MELTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  val.m = melt(val_nmat,  id.vars=row.names)
  new.m = melt(new_nmat, id.vars=row.names)
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##cohort subset and thresholds of percentiles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  comb_cohorts = c()
  samp_cohorts = c()
  for (c in seq(1:length(cohs))){
    temp_coh = val.m[val.m$Var2 %in% make.names(cohs[[c]]),]
    print(length(table(as.character(temp_coh$Var2))))
    temp_coh$ab = names(cohs[c])
    comb_cohorts = rbind(comb_cohorts, temp_coh)
  }
  comb_cohorts$ab_Var1 = paste0(comb_cohorts$ab, "_", comb_cohorts$Var1)
  samp_cohorts = new.m
  ##Get percentiles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  
  #~ ANTIBODY THRESHOLDING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # if signal for antibody below sample igg then turned to minimum of validation  cohort
  # if signal above do nothing
  samp_cohorts['detectable']=TRUE
  newsamp = comb[,paste0("newbatch", "__", (samp)),drop=F]
  rbigg = newsamp[which(rownames(newsamp)=="RbAb-IgG"),]
  mmigg = newsamp[which(rownames(newsamp)=="MmAb-IgG1"),]
  
  for (i in seq(1:length(ab_order))){
    ab = ab_order[i]
    idx = which(samp_cohorts$Var1==ab)
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
      samp_cohorts[idx,"detectable"] = FALSE
      samp_cohorts[idx,"newvalue"] = minval
    }
    else{
      samp_cohorts[idx,"newvalue"] = samp_cohorts[idx,"value"]
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ BINNING FOR HEATMAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  thresh=data.frame("5_95"=c("5%", "95%"), "15_85"=c("15%", "85%"), "25_75"=c("25%","75%"), stringsAsFactors = F)
  cohorder = c()
  for (cidx in seq(1:length(coh_percs))){
    print(names(coh_percs)[cidx])
    coh = names(coh_percs)[cidx]
    for (ab in levels(comb_cohorts$Var1)){
      idxab = which(samp_cohorts$Var1==ab)
      cohtemp = coh_percs[[coh]][ab,]
      for (tidx in seq(1:ncol(thresh))){
        lowthresh = cohtemp[thresh[1,tidx]]
        highthresh = cohtemp[thresh[2,tidx]]
        classname = paste0(names(coh_percs)[cidx], "__", colnames(thresh)[tidx])
        samp_cohorts[idxab,classname] = ifelse(samp_cohorts[idxab,"detectable",drop=T]==FALSE, yes = "ND", 
                                               no = ifelse(samp_cohorts[idxab,"value",drop=T] <= lowthresh, yes="low",
                                                           no = ifelse(samp_cohorts[idxab,"value",drop=T] >= highthresh, yes="high", 
                                                                       no="indeterminate")))
        if(!classname %in% cohorder) {cohorder = c(cohorder, classname)}
      }
    }
  }
  
    #PLOTTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    msamp = melt(samp_cohorts, id.vars = list("Var1", "Var2"), measure.vars = which(colnames(samp_cohorts) %in% cohorder))
    msamp$sampnames = paste0(make.names(sapply(strsplit(as.character(msamp$Var2), "__"), `[`, 2)))
    msamp$cohortnames = paste0(make.names(sapply(strsplit(as.character(msamp$variable), "__"), `[`, 1)))
    msamp$names = paste0(msamp$cohortnames,"__", msamp$sampnames)
    
    #~ binmap
    hmlist <- list()
    for (thresholds in colnames(thresh)){
      mmsamp = msamp[grepl(thresholds, msamp$variable),]
      mmsamp$Var1 = factor(mmsamp$Var1, levels=ab_order)
      mmsamp$cohortnames = factor(mmsamp$cohortnames, levels = unique(mmsamp$cohortnames))
      #mmsamp$Var2 = factor(mmsamp$Var2, levels=cohorder)
      hm= ggplot(mmsamp) +
        geom_tile(aes(x=cohortnames, y=Var1, fill=factor(value)),colour = "grey50") +
        scale_fill_manual(values=c("indeterminate"="white", "ND"="black","low"="yellow","high"="blue")) +
        scale_x_discrete(name ="Cohort", breaks = mmsamp$cohortnames, labels = mmsamp$names) +
        labs(x="Cohort", title=paste0("Threshold on ", thresh[1,thresholds], " and ", thresh[2,thresholds],"\n", unique(mmsamp$sampnames)), y="Antibody") +
        theme(axis.text.x = element_text(size=6, colour=c("black"), angle = 90, hjust=1, vjust=0.5),
              axis.text.y = element_text(size=6),
              axis.title.x = element_text(size=6),
              axis.title.y = element_text(size=6),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray60"),
              panel.background = element_blank(),
              plot.title = element_text(hjust = 0.5, size=9)) 
        
      hmlist[[thresholds]] = hm
    } 
    
    pdf(file=sprintf("%s_binmaps.pdf", samp), width=12, height=7)
    grid.arrange(hmlist[["X5_95"]], hmlist[["X15_85"]], hmlist[["X25_75"]], ncol=3)
    dev.off()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Boxplot
    
    samp_cohorts_box = c()
    for (coh in unique(comb_cohorts$ab)){
      samp_cohorts_box = rbind(samp_cohorts_box, data.frame(samp_cohorts, "ab"=coh))}
    
    #~ BOXPLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    palette = c("#FF0000FF","#004CFFFF","#00A600FF","#984ea3","#ff7f00","#a65628")
    
    pwboxlist = list()
    
    for (pathway in as.character(unique(pathways$Pathway))){
      if (!pathway %in% c("Non-specific")){
        antibodies = as.character(pathways$ab_ref.X.AbID[pathways$Pathway==pathway])
        comb_cohorts_box = comb_cohorts[comb_cohorts$Var1 %in% antibodies,]
        samp_cohorts_boxab = samp_cohorts_box[samp_cohorts_box$Var1 %in% antibodies,]
        
        bp = ggplot(comb_cohorts_box, aes(x=ab, y=value)) + 
          geom_boxplot() +
          facet_wrap(~Var1, scale="free", nrow=1) +
          geom_point(data=samp_cohorts_boxab[samp_cohorts_boxab$detectable==TRUE,], mapping=aes(x=ab, y=newvalue, colour=Var2), shape=8, size=2) +
          geom_text(data=samp_cohorts_boxab[samp_cohorts_boxab$detectable==FALSE,], mapping=aes(x=ab, y=newvalue, colour=Var2), label="ND", size=3, position=position_jitter(width=c(0.03))) +
          labs(y="log batch corrected counts", title=pathway) +
          scale_color_manual(values=palette[1:length(samp)]) +
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
    lay <- rbind(c(1,1,1,1,NA,NA,NA,NA),
                 c(2),
                 c(3),
                 c(4,4,4,NA,5,5,NA,NA))
    
    
    pdf(file=sprintf("%s_boxplots.pdf", samp), width=10, height=8.5)
    
    grid.arrange(grobs=pwboxlist, layout_matrix=lay, top=samp)
    dev.off()
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    write.table(lcomb[!grepl("NEG|POS", colnames(lcomb)),], file=paste0("ruv_figures/",samp,"_lograw.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
    write.table(t(RUVcorrected), file=paste0("ruv_figures/",samp,"_ruvcorrected.txt"),sep="\t", quote = F, row.names = T, col.names = NA)
  
  }
  
tar(tarfile=paste0("ruv_figures.tar.gz"), files=paste0("ruv_figures"), compression="gzip", tar="tar")

  
  
  
  
  
  

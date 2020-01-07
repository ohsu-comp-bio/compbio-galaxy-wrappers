#!/usr/bin/env Rscript

##
## QC metrics, RUV and comparison with MBC cohort 
## Usage
## ./ruv_mbc.R -i 3_IGG_NORMALIZED.tsv --validation_file validation_samples_normalized.txt --mbc_md_file validation_mbc_metadata.txt --ab_ref_file ANTIBODY_REFERENCE.csv
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ruv))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gridExtra))


## ARGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
version="2.0"


parser <- ArgumentParser()
parser$add_argument("-i", help="3_IGG_NORMALIZED.tsv", dest="input_file")
parser$add_argument("--validation_file", type="character", default="/Users/patterja/Box\ Sync/NANOSTRING/REFERENCE_FILES/validation_samples_sub_normalized.txt",
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--mbc_md_file", type="character", default= "/Users/patterja/Box Sync/NANOSTRING/REFERENCE_FILES/validation_mbc_metadata.txt",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--ab_ref_file", type="character", default= "/Users/patterja/Box\ Sync/NANOSTRING/REFERENCE_FILES/ANTIBODY_REFERENCE.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--include_bad_Ab", action="store_true", default=FALSE,
                    dest="include_bad_Ab", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
includeBCCL = args$includeBCCL

input_file = args$input_file
validation_file = args$validation_file
mbc_md_file = args$mbc_md_file
ab_ref_file = args$ab_ref_file
include_bad_Ab = args$include_bad_Ab


## CONTROLS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF")
ctrlregex = gsub("\\+", "\\\\+", paste0(paste0(controls, collapse = "|"),"|", "MDA468"))

## FUNCTIONS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pcaplt <- function (mat, title = "PCA Plot", repdf) {
  var = mat[apply(mat, 1, var, na.rm = TRUE) != 0, ]
  cc.var = var[complete.cases(var), ]
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  PC1_and_PC2 = data.frame(PC1 = pca_prcomp$x[, 1], PC2 = pca_prcomp$x[,2], type = rownames(pca_prcomp$x))
  perc = (pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) * 100
  labs <- sapply(seq_along(perc), function(i) {
    paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  PCsmd = cbind(PC1_and_PC2, repdf[match(PC1_and_PC2$type, repdf$sampcolumn), 
                                   c("samp", "batch")])
  p = ggplot(PCsmd,aes_string("PC1", "PC2", col="batch")) + 
    geom_point(size = 1.5) + 
    geom_text(aes(label = type), vjust = -1, size=2) + 
    labs(title = title,x = labs[1], y = labs[2]) + 
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "gray"), 
          plot.title = element_text(hjust = 0.5), 
          legend.text = element_text(size = 4), legend.position = "right")
  return(p)
}
## OUTPUT PREP ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#dir.create(paste0(dirname(normalizePath(input_file)), "/OUTPUT"), showWarnings = F)
#output_dir = paste0(dirname(normalizePath(input_file)), "/OUTPUT")
output_dir = "."
#dir.create(file.path(paste0(dirname(normalizePath(input_file)), "/OUTPUT/ruv_figures")), showWarnings = F)
dir.create("ruv_figures", showWarnings = F)

## INPUT DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (include_bad_Ab){
  print("keeping all antibodies: IgG and antibodies that did not perform well or did not have dynamic range")
  #validation file
  igg_geosamp = read.csv(validation_file, sep= "\t", row.names = 1, check.names = F)
  colnames(igg_geosamp)=make.names(gsub("\\.1$", "", colnames(igg_geosamp)))
  igg_geosamp=igg_geosamp[,grepl("062519", colnames(igg_geosamp))]
  ligg_geosamp=log2(igg_geosamp+1)
  #AB FILE
  ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)
  #NEW BATCH
  new_batch = data.matrix(read.table(file = input_file, sep="\t", row.names=1, stringsAsFactors=F, header=T))
  #MBC 
  mbc_md= read.table(file = mbc_md_file, sep="\t", row.names=1, stringsAsFactors=F, header=T)
  
} else {
  omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2")
  omitregex = paste0(paste0("^", omit), collapse = "|")
  #validation file
  igg_geosamp = read.csv(validation_file, sep= "\t", row.names = 1, check.names = F)
  igg_geosamp=igg_geosamp[,!grepl("062519", colnames(igg_geosamp))]
  igg_geosamp = igg_geosamp[!grepl(omitregex, rownames(igg_geosamp)),]
  
  
  ligg_geosamp=log2(igg_geosamp+1)
  colnames(ligg_geosamp)=make.names(gsub("\\.1$", "", colnames(ligg_geosamp)))

  #AB FILE
  ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)
  ab_ref = ab_ref[!grepl(omitregex, ab_ref$X.AbID),]
  #NEW BATCH
  new_batch = data.matrix(read.table(file = input_file, sep="\t", row.names=1, stringsAsFactors=F, header=T))
  new_batch = new_batch[!grepl(omitregex, rownames(new_batch)),]
  #MBC 
  mbc_md= read.table(file = mbc_md_file, sep="\t", row.names=1, stringsAsFactors=F, header=T)
}

#AB_ORDER
ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
ab_order = ab_order[!grepl(omitregex, ab_order)]

#VALIDATION QC MATRIX #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lctl_norm = as.matrix(ligg_geosamp)[,grepl(ctrlregex, colnames(ligg_geosamp))]
lctl_norm.m = melt(lctl_norm)
lctl_norm.m$batch = sapply(strsplit(as.character(lctl_norm.m$Var2), "__"), `[`, 1)
lctl_norm.m$ctlname = sapply(strsplit(as.character(lctl_norm.m$Var2), "__"), `[`, 2)
lctl_norm.m$ctl_probe = factor(paste0(lctl_norm.m$ctlname,"_", lctl_norm.m$Var1))


validation_stats = data.frame(
  #do.call(rbind, (tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, summary))),  #quartiles
  "min" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, min),
  "mean" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean),
  "max" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, max),
  "stddev" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd),
  "coeff_var" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)/tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean),
  "mean.minus.2sd" =tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean) - 
    2*(tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)),
  "mean.plus.2sd" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean) + 
    2*(tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)))

#only use controls 
ctl_newbatch = melt(log2(new_batch[,make.names(controls)]+1))
ctl_newbatch$ctl_probe = paste0(gsub("\\.1$", "", sapply(strsplit(as.character(ctl_newbatch$Var2), split = "_"), tail,1)),"_", ctl_newbatch$Var1)
ctl_newbatch$cellline = factor(gsub("\\.1$", "", sapply(strsplit(as.character(ctl_newbatch$Var2), split = "_"), tail,1)))

ctl_validation =cbind(validation_stats, "new_batch"=ctl_newbatch[match(rownames(validation_stats),ctl_newbatch$ctl_probe), "value"])
ctl_validation$status = ifelse(ctl_validation$new_batch < ctl_validation$mean.minus.2sd, 
                               (ctl_validation$mean - ctl_validation$new_batch)/ctl_validation$stddev, 
                               ifelse(ctl_validation$new_batch >= ctl_validation$mean.plus.2sd, 
                                      (ctl_validation$new_batch - ctl_validation$mean)/ctl_validation$stddev, "PASS"))

write.table(ctl_validation, file = paste0(output_dir, "/", "qc_controls.tsv"), 
            sep="\t", quote = F, row.names = T, col.names = NA)

## QC LINEAR MODEL PLOTS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Comparing newQC with old QC 
#simple linear regression analysis

lmfitqc = lm(new_batch ~ mean, data=ctl_validation)

#~ plot linear model
ggReg = ggplot(lmfitqc$model, aes_string(x = names(lmfitqc$model)[2], y = names(lmfitqc$model)[1])) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Antibody Agreement: point per cell line per antibody",
                     "\nCorrelation= ",signif(cor((lmfitqc$model)[1],(lmfitqc$model)[2]), 5),
                     "\nAdj R2= ",signif(summary(lmfitqc)$adj.r.squared, 5),
                     "\nIntercept=",signif(lmfitqc$coef[[1]],5),
                     " p=",signif(summary(lmfitqc)$coef[1,4], 5),
                     "\nSlope=",signif(lmfitqc$coef[[2]], 5),
                     " p=",signif(summary(lmfitqc)$coef[2,4], 5)), 
       x="mean of validation samples")

twosd = 2*sd(ctl_validation$new_batch-ctl_validation$mean)
outdat =  ctl_validation[which((ctl_validation$new_batch-ctl_validation$mean) < -twosd|(ctl_validation$new_batch-ctl_validation$mean) > twosd),]
outdat$cellline = sapply(strsplit(as.character(rownames(outdat)), "_"), `[`, 1)
outdat$probe = sapply(strsplit(as.character(rownames(outdat)), "_"), `[`, 2)

baplot_ctl = ggplot(ctl_validation) +
  geom_point(aes(x=mean, y=new_batch-mean)) +
  geom_hline(yintercept=c(0, twosd, -twosd), color="red", linetype = 2) +
  geom_text(data=data.frame(x=0,y=c(-twosd, twosd)), aes(x, y), label = c("2SD", "-2SD"), color="red") +
  geom_point(data = outdat, aes(x=mean,  y=outdat$new_batch-outdat$mean,color=outdat$cellline)) +
  #geom_text(data = outdat, aes(x=mean,  y=outdat$new_batch-outdat$mean,
  #                             label=rownames(outdat), hjust=-0.1), size=3) +
  labs(title = "Mean-Difference Plot: labeled if +/- >2SD", x="Mean of validation batches", y="New Batch - Mean of validation batches")

baplot_ab = ggplot(ctl_validation) +
  geom_point(aes(x=mean, y=new_batch-mean)) +
  geom_hline(yintercept=c(0, twosd, -twosd), color="red", linetype = 2) +
  geom_text(data=data.frame(x=0,y=c(-twosd, twosd)), aes(x, y), label = c("2SD", "-2SD"), color="red") +
  geom_point(data = outdat, aes(x=mean,  y=outdat$new_batch-outdat$mean,color=outdat$probe)) +
  labs(title = "Mean-Difference Plot: labeled if +/- >2SD", x="Mean of validation batches", y="New Batch - Mean of validation batches")


print("saving QC plot")
ggsave(filename =paste0(output_dir,"/QC_antibody_linear_plot.pdf"), device = "pdf", plot=ggReg,width = 8, height = 8)
ggsave(filename =paste0(output_dir,"/QC_antibody_meandiff_plot.pdf"), device = "pdf", grid.arrange(baplot_ab, baplot_ctl),width = 8, height = 5)


## RUV DATASET ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#for sample in newbatch
mbc_percentile=data.frame(row.names = rownames(ligg_geosamp))

for (samp in setdiff(colnames(new_batch), make.names(controls))) {
  print(samp)
  #COMBINING: rows combined with mbcctl, columns=ctrls and 1 samp only
  newsampctl = new_batch[,c(grep(ctrlregex, colnames(new_batch)),which(colnames(new_batch)==samp))]
  colnames(newsampctl) = paste0("newbatch", "__", colnames(newsampctl))
  combined_norm = cbind(igg_geosamp, newsampctl[match(rownames(igg_geosamp), rownames(newsampctl)),])
  combined_md = data.frame(batch = sapply(strsplit(as.character(colnames(combined_norm)), "__"), `[`, 1), 
                           samp = gsub("\\.1$", "", sapply(strsplit(as.character(colnames(combined_norm)), "__"), `[`, 2)),
                           sampcolumn = c(colnames(combined_norm)),stringsAsFactors = F)
  combined_md$mbc= ifelse(combined_md$batch %in% mbc_md$BatchID, #MBC batch only
                          yes=ifelse(combined_md$samp %in% make.names(controls), yes="MBCcontrol", #controls MBC 
                                     no=ifelse(combined_md$samp %in% make.names(mbc_md$Sample), yes="MBCsample",no="notMBC")),#in mbc md sample
                          no="notMBC")
  #rep matrix based only on controls in MBC batches
  combined_md$reps=ifelse(combined_md$mbc=="MBCcontrol", yes=combined_md$samp, no=combined_md$sampcolumn)
  #~ replicated matrix 
  repmat = replicate.matrix(combined_md$reps)
  
  #~ plotting pre-ruv figures
  lctls = log2(as.matrix(combined_norm)[,grep(ctrlregex, colnames(combined_norm))]+1)
  
  limits_rle = max(as.matrix(lctls))-median(as.matrix(lctls))
  rle_orig = ruv_rle(Y = t(lctls), rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% colnames(lctls),]), 
                ylim=c(-limits_rle,limits_rle)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=samp), angle=90, hjust=0, size=2)+
    theme(legend.position = "right", legend.text = element_text(size=6)) + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression \nNormalized Data (IGG corrected) \n",  samp))
  pca_orig = pcaplt(lctls, title = paste0("log2 (Normalized Signal + 1)\n", samp), repdf=combined_md)
  
  #~ RUV
  RUVcorrected = RUVIII(Y=t(log2(combined_norm +1)), ctl=c(1:nrow(combined_norm)), M=repmat, k=1)
  #~ plotting RUV processed figures
  ctl_ruv = t(RUVcorrected)[,grep(ctrlregex, colnames(t(RUVcorrected)))]
  rle_ruv = ruv_rle(Y = RUVcorrected[grep(ctrlregex, rownames(RUVcorrected)),], 
                    rowinfo = as.matrix(combined_md[combined_md$sampcolumn %in% colnames(ctl_ruv),]), 
                    ylim=c(-limits_rle,limits_rle)) +
    geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
    geom_text(aes(x = rle.x.factor, y=-limits_rle-0.5, label=samp), angle=90, hjust=0, size=2) +
    theme(legend.position = "right") + 
    labs(colour = "batch") +
    geom_hline(yintercept = 0, 
               linetype = "dotted", colour = "cyan") + 
    ggtitle(paste0("Relative Log Expression RUV processed\n", samp))
  pca_ruv = pcaplt(ctl_ruv, title = paste0("RUV Processed\n", samp), repdf=combined_md)
  
  #~ split these apart makes plotting easier
  #mbc only and newsamp only
  mbc_ruv = t(RUVcorrected)[,as.character(combined_md$sampcolumn[combined_md$mbc=="MBCsample"])]
  new_ruv = data.frame("Var1" = rownames(t(RUVcorrected)),
                       "Var2"= samp,
                       "value" = t(RUVcorrected)[rownames(t(RUVcorrected)), paste0("newbatch__", samp)])
  
  #~ Get percentile using ecdf
  #~ melt for plotting
  mbcruv.m = melt(mbc_ruv,  id.vars=row.names)
  mbc_ecdf = tapply(mbcruv.m$value, mbcruv.m$Var1, ecdf)
  mbc_percentile[,samp] =  as.vector(apply(new_ruv, 1, function(x) 
    ((mbc_ecdf[[as.character(x["Var1"])]](x[["value"]]))))[rownames(mbc_percentile)])
  
  
  #max and min
  ruv = as.matrix(t(RUVcorrected))
  mruv = melt(ruv)
  
  ruv_stats = data.frame(
    #do.call(rbind, (tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, summary))),  #quartiles
    "min" = tapply(mruv$value, mruv$Var1, min),
    "q1" = tapply(mruv$value, mruv$Var1, function(x) quantile(x, 0.25)),
    "mean" = tapply(mruv$value, mruv$Var1, median),
    "q3" = tapply(mruv$value, mruv$Var1, function(x) quantile(x, 0.75)),
    "max" = tapply(mruv$value, mruv$Var1, max))
  
  
  mbcruv.m$Var1 = factor(mbcruv.m$Var1, levels=ab_order)
  new_ruv$Var1 = factor(new_ruv$Var1, levels=ab_order)
  ruv_stats$ab = factor(rownames(ruv_stats), levels=ab_order)
  
  
  mbc_labels= paste0(as.character(levels(new_ruv$Var1))," (",round(new_ruv$value, 1), ",",round(mbc_percentile[as.character(levels(new_ruv$Var1)), samp]*100,0),")")
 
  #~ boxplot of MBC

  bp_mbcruv = ggplot(data.frame(mbcruv.m)) +
    geom_point(data=new_ruv, mapping=aes(x=factor(Var1), y=value), colour=c('red'), shape=8, size=2) +
    geom_linerange(
      data=ruv_stats, aes(x=ab, ymin = min, ymax = max),
      color = "#808080", 
      size = 7, 
      alpha = 0.7) +
    #geom_point(data = ruv_stats, aes(x=ab, y=max),shape=93, fill="grey") +
    geom_boxplot(aes(factor(Var1), as.numeric(value)), outlier.colour = NA) +
    geom_point(data=new_ruv, mapping=aes(x=factor(Var1), y=value), colour=c('red'), shape=8, size=2) +
    scale_x_discrete(labels=paste0(as.character(levels(new_ruv$Var1))," (",
                                   round(mbc_percentile[as.character(levels(new_ruv$Var1)), samp]*100,0),")")) +
    labs(x="Antibody (Percentiles)", title=paste0(samp,  "\n within Distribution of Metastatic Breast Cancers"), y="RUVnormalized") +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5, vjust=0),
          legend.text=element_text(size=8),
          legend.position="right",
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black")) +
    coord_flip()
  
  print("saving RUV plots")
  #~ save all plots
  ggsave(file=paste0(output_dir, "/",samp, "_MBC.pdf"), bp_mbcruv, device="pdf", width = 8, height = 6)
  ggsave(file=paste0(output_dir, "/ruv_figures/",samp, "_controls_norm_RLE.pdf"), 
         device="pdf", rle_orig, width = 9, height = 4.5)
  ggsave(file=paste0(output_dir, "/ruv_figures/",samp, "_controls_norm_PCA.pdf"), 
         device="pdf", pca_orig, width = 8, height = 7)
  ggsave(file=paste0(output_dir, "/ruv_figures/",samp, "_controls_RUV_RLE.pdf"), 
         device="pdf", rle_ruv, width = 9, height = 4.5)
  ggsave(file=paste0(output_dir, "/ruv_figures/",samp, "_controls_RUV_PCA.pdf"), 
         device="pdf", pca_ruv, width = 8, height = 7)
  
  
  write.table(x = t(RUVcorrected), file=paste0(output_dir, "/ruv_figures/",samp, "_RUVcorrected.tsv"), 
              sep="\t", quote = F, row.names = T, col.names = NA)
}

tar(tarfile=paste0(output_dir,"/ruv_figures.tar.gz"), files=paste0(output_dir, "/ruv_figures"), compression="gzip", tar="tar")
write.table(x = round(mbc_percentile, 2), file=paste0(output_dir, "/", "MBC_percentiles.tsv"), sep="\t", quote = F, row.names = T, col.names = NA)


  
  
  
  

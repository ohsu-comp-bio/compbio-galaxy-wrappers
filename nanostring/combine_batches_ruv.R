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
suppressPackageStartupMessages(library(nanoprot))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))


version="4.0"

parser <- ArgumentParser()
parser$add_argument("--data_dir", type="character", default="/Volumes/OHSU/CLINICAL/Nanostring/output", 
                    dest="data_dir", help="directory name")
parser$add_argument("--comb_sheet", type="character", dest="comb_sheet", help="ruv samplesheet")
parser$add_argument("--md_file", type="character", default= "/Users/patterja/Box/NANOSTRING/nanostring_metadata.xlsx",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--ab_ref_file", type="character", default= "/Volumes/OHSU/CLINICAL/Nanostring/REFERENCE_FILES/ANTIBODY_REFERENCE_v1.0.csv",
                    dest="ab_ref_file", help="ANTIBODY_REFERENCE.csv")
parser$add_argument("--include_ctrls", action="store_true", default=FALSE,
                    dest="include_ctrls", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
data_dir = args$data_dir
comb_sheet = args$comb_sheet
md_file = args$mbc_md_file
ab_ref_file = args$ab_ref_file
include_ctrls = args$include_ctrls

#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468.control","MDA468+EGF")
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2")
omitregex = paste0(paste0("^", omit), collapse = "|")
ab.ctrl = "IgG|POS|NEG|^S6|^Histone"



#FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#outdirectory. Not sure how this will integrate with galaxy
outdir = dirname(normalizePath(comb_sheet))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~ samplesheet
print(paste0("output files will be in current dir: ", outdir))

samps2batchcorr = read.table(file= comb_sheet, sep="\t", stringsAsFactors = F, header=T)
sampname = gsub("_samplesheet.txt", "", basename(comb_sheet))

#~ metadata
md = read.xlsx(xlsxFile=md_file, sheet="nanostring_metadata", check.names=T)
md$sampcolumn = make.names(paste0(md$Batch, "__", md$Sample.Name))

# antibody metadata 
ab_ref = read.csv(ab_ref_file, sep=",", stringsAsFactors=F)
pathways = data.frame(ab_ref$X.AbID, Pathway = sapply(strsplit(as.character(ab_ref$Pathway), ","), `[`, 1))
pathways = rbind(pathways, data.frame(ab_ref$X.AbID, Pathway=sapply(strsplit(as.character(ab_ref$Pathway), ","), `[`, 2)))
pathways = rbind(pathways, data.frame(ab_ref$X.AbID, Pathway=sapply(strsplit(as.character(ab_ref$Pathway), ","), `[`, 3)))

pathways = pathways[complete.cases(pathways),]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ PARSE ALL SAMPLES FROM ALL BATCHES IN SAMPSHEET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
raw_dats = c()
raw_samps = c()
for (f in 1:length(unique(samps2batchcorr$Batch))){
  
  #get raw data
  batch_name = unique(samps2batchcorr$Batch)[f]
  file_name = as.character(paste0(data_dir,"/",batch_name ,"/rawdata.txt"))
  if (file.exists(file_name)){
    batch = read.table(file = file_name, sep="\t", row.names=2, stringsAsFactors=F, header=T, check.names = T)
  } else {
    
    stop(sprintf("Batch QC failed. No rawdata.txt file exists! Batch %s", batch_name))
  }
  batch[,c("CodeClass", "Accession")] <- NULL
  
  ## TEST ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #TEST to see if sample names in batch directories
  ss_samps=make.names(samps2batchcorr[which(samps2batchcorr$Batch==batch_name),"Sample.Name"])
  for (samp in ss_samps){
    if(samp %in% colnames(batch)){
      print(paste0("Sample:", samp, " in ", file_name))
    } else {
      print(paste0("Sample:", samp, " NOT in ", file_name, "\n Check sample name."))
      stop()
    }
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
  
  #get only samples you want to batch correct and the controls
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


# COMBINING DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## combine for replicate structure: rows combined with validation, ctrls and samples
lcombined = log(raw_dats+1)

# COMBINED DATA METADATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create a new metadata for replicate structure
combined_md = raw_samps
combined_md$deid = ifelse(combined_md$sample %in% make.names(controls), yes=combined_md$sample, no="sample")



#~ RUV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
combined_md$reps = ifelse(combined_md$sample %in% make.names(controls), yes=combined_md$sample, no=combined_md$fullname)

repmat = replicate.matrix(combined_md$reps)
idx_ctl = grep(ab.ctrl, rownames(lcombined))

RUVcorrected = RUVIII(Y=t(lcombined), ctl=idx_ctl, M=repmat, k=2, include.intercept = FALSE, inputcheck = F)
bcdat = t(RUVcorrected[,!grepl("POS|NEG", rownames(lcombined))])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## because i got sick of calling it all the time
samps = raw_samps$fullname[!raw_samps$sample %in% make.names(controls)]
celllines = setdiff(raw_samps$fullname, samps)


#~ QC PLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("qc_figures", showWarnings = F)

#~ RLE box plot
lctls=lcombined[!grepl("POS|NEG", rownames(lcombined)),celllines]

limits_rle_raw = max(as.matrix(lctls))-median(as.matrix(lctls))
rle_rawbox = ruv_rle(Y = t(lctls), rowinfo = as.matrix(combined_md[combined_md$fullname %in% celllines,]), 
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
                   rowinfo = as.matrix(combined_md[combined_md$fullname %in% celllines,]), 
                   ylim=c(-limits_rle_bc,limits_rle_bc)) +
  geom_point(aes(x = rle.x.factor, y = middle, colour = batch)) +
  geom_text(aes(x = rle.x.factor, y=-limits_rle_bc-0.5, label="RUVcorrected"), angle=90, hjust=0, size=2) +
  theme(legend.position = "right") + 
  labs(colour = "batch") +
  geom_hline(yintercept = 0, 
             linetype = "dotted", colour = "cyan") + 
  ggtitle(paste0("Relative Log Expression RUV processed\n", sampname, "_", "RUVcorrected"))



# NORMAL HEATMAP
pdf(paste0("qc_figures/",sampname,"_signal_heatmap.pdf"), width = 7, height = 8)
## rle raw
rownames(combined_md) = combined_md$fullname
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
pdf(paste0("qc_figures/",sampname,"_RLE_heatmap.pdf"), width = 7, height = 8)
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
pcacol= combined_md$deid[match(colnames(lctls), combined_md$fullname)]
pca_raw = pcaplot(mat = lctls, title = paste0("log2 (Normalized Signal + 1)\n", sampname, "_", "PRE-RUV"), col=pcacol)

pcacol= combined_md$deid[match(colnames(ctl_bc), combined_md$fullname)]
pca_ruv = pcaplot(ctl_bc, title = paste0("RUV Processed\n", sampname, "_", "RUVcorrected"), pcacol)

#~ TRA 
#t(RUVcorrected) is post norm, comb is pre norm, using no antibody filtered data
eruv = t(exp(ctl_bc))
rawctls = t(raw_dats)[rownames(eruv),colnames(eruv)]
tra = matrix(NA, nrow = 0, ncol = 4)

for (ctrl in make.names(controls)){
  #get matrix of celllines for controls not including batches of interest
  valid.ctrl_names = combined_md$fullname[combined_md$samp == ctrl]

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

allselraw = t(raw_dats[!grepl("POS|NEG", rownames(lcombined)),celllines])
allselruv = t(exp(bcdat[,celllines]))

for (ctrl in make.names(controls)){
  #get matrix of celllines for controls not including batches of interest
  valid.ctrl_names = celllines[grep(ctrl, celllines)]
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
  labs(title=paste0("Distribution of TRA between RAW and RUV Corrected between Controls of Batches being compared\n"),
       y="log(validation count/new batch count)") +
  theme(legend.position = "bottom")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


print("saving batch correction plots")
pdf(file=paste0("qc_figures/",sampname,"_TRA.pdf"), width = 8, height = 8)
print(ptra)
print(ptra_btwn)
dev.off()

pdf(paste0("qc_figures/",sampname,"_RLE_boxplots.pdf"), width = 8, height = 4)
print(rle_rawbox)
print(rle_ruvbox)
dev.off()

pdf(paste0("qc_figures/",sampname,"_PCA.pdf"), width = 8, height = 8)
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
  pathways = pathways[pathways$ab_ref.X.AbID %in% ab_order,]
} else {
  fbcdat = bcdat[!grepl(omitregex, rownames(bcdat)),]
  #AB_ORDER
  ab_order = ab_ref$X.AbID[order(ab_ref$Target)]
  ab_order = ab_order[!grepl(omitregex, ab_order)]
  pathways = pathways[pathways$ab_ref.X.AbID %in% ab_order,]
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ SPLITTING AND MELTING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samps_dat = fbcdat[,!colnames(fbcdat) %in% celllines,drop=F]

samps.datm = data.frame(melt(samps_dat, id.vars=row.names))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ANTIBODY THRESHOLDING setting with detectable flag ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## if signal for antibody below sample igg then turned to minimum of validation  cohort
## if signal above do nothing
samps.datm['detectable']=TRUE
for (samp in samps){
  newsamp = raw_dats[,samp, drop=F]
  #newsamp = lcombined[,samp, drop=F]
  rbigg = newsamp[which(rownames(newsamp)=="RbAb-IgG"),]
  mmigg = newsamp[which(rownames(newsamp)=="MmAb-IgG1"),]
  for (i in seq(1:length(ab_order))){
    ab = (ab_order)[i]
    idx = which(samps.datm$Var1==ab & samps.datm$Var2==samp)

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
      samps.datm[idx,"detectable"] = FALSE
      samps.datm[idx,"newvalue"] = NA
    } else{
      samps.datm[idx,"newvalue"] = samps.datm[idx,"value"]
    }
  }
}

#~ ORDERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samp_box = samps.datm
samp_box$Var1 = factor(samp_box$Var1, levels=ab_order)

#~ BOXPLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~ HEATMAP #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samps.datc = acast(samps.datm, Var1 ~ Var2, value.var = "newvalue")
samps.datc = samps.datc[ab_order,sort(colnames(samps.datc))]
samps.datraw = lcombined[!grepl("POS|NEG", rownames(lcombined)),!colnames(lcombined) %in% celllines,drop=F]
samps.datraw = samps.datraw[ab_order,sort(colnames(samps.datraw))]
pdf(file=paste0(sampname, "_comparison_heatmap.pdf"), width = 7, height = 11)

pheatmap(mat = samps.datraw,
         color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
         cluster_cols      = F,
         cluster_rows      = F,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         na_col            = c("white"),
         border_color      = c("black"),
         main              = sprintf("log raw: %s", sampname))


pheatmap(mat = samps.datc,
         color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
         cluster_cols      = F,
         cluster_rows      = F,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         na_col            = c("white"),
         border_color      = c("black"),
         main              = sprintf("RUV corrected: %s", sampname))

samps.datc.na = samps.datc
samps.datc.na[is.na(samps.datc.na)] =-1
pheatmap(mat = samps.datc.na,
         color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
         cluster_cols      = F,
         cluster_rows      = T,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         na_col            = c("white"),
         border_color      = c("black"),
         main              = sprintf("RUV corrected: %s \n Clustered Antibodies", sampname))

pheatmap(mat = samps.datc.na,
         color             = colorRampPalette( c("green", "black", "red"), space="rgb")(10),
         cluster_cols      = T,
         cluster_rows      = T,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         na_col            = c("white"),
         border_color      = c("black"),
         main              = sprintf("RUV corrected: %s \n Clustered Antibodies and Clustered Samples", sampname))

dev.off()

pdf(file=paste0(sampname, "_comparison_heatmap2.pdf"), width = 7, height = 11)

pheatmap(mat = samps.datraw,
         color             = c(colorRampPalette(brewer.pal(9,"YlOrRd"))(9),"#67001f"),
         cluster_cols      = F,
         cluster_rows      = F,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         na_col            = c("white"),
         border_color      = c("black"),
         main              = sprintf("log raw: %s", sampname))

pheatmap(mat = samps.datc,
         color             = c(colorRampPalette(brewer.pal(9,"YlOrRd"))(9),"#67001f"),
         cluster_cols      = F,
         cluster_rows      = F,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         na_col            = c("white"),
         border_color      = c("black"),
         main              = sprintf("RUV corrected: %s", sampname))

samps.datc.na = samps.datc
samps.datc.na[is.na(samps.datc.na)] =-1
pheatmap(mat = samps.datc.na,
         color             = c("#FFFFFF", colorRampPalette(brewer.pal(9,"YlOrRd"))(9),"#67001f"),
         cluster_cols      = F,
         cluster_rows      = T,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         na_col            = c("white"),
         border_color      = c("black"),
         main              = sprintf("RUV corrected: %s \n Clustered Antibodies", sampname))

pheatmap(mat = samps.datc.na,
         color             = c("#FFFFFF", colorRampPalette(brewer.pal(9,"YlOrRd"))(9),"#67001f"),
         cluster_cols      = T,
         cluster_rows      = T,
         show_colnames     = TRUE,
         show_rownames     = TRUE,
         fontsize          = 11,
         na_col            = c("white"),
         border_color      = c("black"),
         main              = sprintf("RUV corrected: %s \n Clustered Antibodies and Clustered Samples", sampname))

dev.off()

#~ TABLES #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(x = t(bcdat), file=paste0(sampname, "_RUVcorrected.csv"),sep=",", quote = F, row.names = T, col.names = NA)
write.table(x = t(raw_dats), file=paste0(sampname, "_raw.csv"), sep=",", quote = F, row.names = T, col.names = NA)


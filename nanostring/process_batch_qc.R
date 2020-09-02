#!/usr/bin/env Rscript

## QC metrics and tables
## Usage: see combined_batches_ruv.R --help
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(jsonlite))

## ARGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
version="2.0"
matrix_version = "20200512"

parser <- ArgumentParser()
parser$add_argument("-i", help="rawdata.txt", dest="input_file")
parser$add_argument("--validation_file", type="character", 
                    default=paste0("/Volumes/OHSU/CLINICAL/Nanostring/REFERENCE_FILES/validation_samples_rawdata_", matrix_version,".txt"),
                    dest="validation_file", help="validation file for controls comparison")
parser$add_argument("--pos_file", type="character", 
                    default=paste0("/Volumes/OHSU/CLINICAL/Nanostring/REFERENCE_FILES/knownpositives_v1.0.txt"),
                    dest="pos_file", help="validation file known positive cellline and antibodies")
parser$add_argument("--md_file", type="character", default= "/Users/patterja/Box/NANOSTRING/nanostring_metadata.xlsx",
                    dest="md_file", help="metadata file")
parser$add_argument("--include_bad_Ab", action="store_true", default=FALSE,
                    dest="include_bad_Ab", help="include all antibodies")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()

input_file = args$input_file
validation_file = args$validation_file
pos_file = args$pos_file
md_file = args$md_file
include_bad_Ab = args$include_bad_Ab


## GLOBAL ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468+EGF")
ctrlregex = gsub("\\+", "\\\\+", paste0(paste0(controls, collapse = "|"),"|", "MDA468"))

## FUNCTIONS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## INPUT DATA ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pos = read.csv(pos_file, sep= "\t", row.names = 3, check.names = F)

#metadata
md = read.xlsx(xlsxFile=md_file, sheet="nanostring_metadata", check.names=T)
md$combined_name = paste0(md$Batch, "__", md$Sample.Name)

#rawdata
new_batch = read.table(file = input_file, sep="\t", row.names=2, stringsAsFactors=F, header=T, check.names = T)
new_batch[,c("CodeClass", "Accession")] <- NULL
if (!any(make.names(controls) %in% colnames(new_batch))){
  stop("no controls in this batch")
}

## log and include_bad_ab flag ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


if (include_bad_Ab){
  print("keeping all antibodies: IgG and antibodies that did not perform well or did not have dynamic range")
  # validation file
  validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = F)
  validation=log2(validation+1)
  
  # NEW BATCH
  
  lnew_batch = log2(new_batch+1)

} else {
  omit = c("p-TSC2", "TSC2", "NEG","POS")
  omitregex = paste0(paste0("^", omit), collapse = "|")
  
  #validation file
  validation = read.csv(validation_file, sep= "\t", row.names = 1, check.names = T)
  validation = validation[!grepl(omitregex, rownames(validation)),]
  lvalidation=log2(validation+1)
 
  #NEW BATCH
  new_batch = new_batch[!grepl(omitregex, rownames(new_batch)),]
  lnew_batch = log2(new_batch+1)
}


#VALIDATION QC MATRIX #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

control_samples = md$combined_name[md$Study %in% c("control")]
lctl_norm = as.matrix(lvalidation[,colnames(lvalidation) %in% make.names(control_samples)])
lctl_norm.m = melt(lctl_norm)
lctl_norm.m$batch = sapply(strsplit(as.character(lctl_norm.m$Var2), "__"), `[`, 1)
lctl_norm.m$ctlname = sapply(strsplit(as.character(lctl_norm.m$Var2), "__"), `[`, 2)
lctl_norm.m$ctl_probe = factor(paste0(lctl_norm.m$ctlname,"__", lctl_norm.m$Var1))


validation_stats = data.frame(
  #do.call(rbind, (tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, summary))),  #quartiles
  "min" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, min),
  "mean" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean),
  "max" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, max),
  "stddev" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd),
  "coeff_var" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)/tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean),
  "mean.minus.2sd" =tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean) - 
    3*(tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)),
  "mean.plus.2sd" = tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, mean) + 
    3*(tapply(lctl_norm.m$value, lctl_norm.m$ctl_probe, sd)))

# NEWBATCH CONTROLS  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#only use controls 
ctl_newbatch = melt(as.matrix(lnew_batch[,make.names(controls)]))
ctl_newbatch$ctl_ab = paste0(ctl_newbatch$Var2,"__", ctl_newbatch$Var1)
ctl_newbatch$cellline = ctl_newbatch$Var2

# COMBINE INTO TABLE  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ctl_validation =cbind(validation_stats, "new_batch"=ctl_newbatch[match(rownames(validation_stats),ctl_newbatch$ctl_ab), "value"])
ctl_validation$status = ifelse(ctl_validation$new_batch < ctl_validation$mean.minus.2sd, 
                               (ctl_validation$new_batch - ctl_validation$mean)/ctl_validation$stddev, 
                               ifelse(ctl_validation$new_batch >= ctl_validation$mean.plus.2sd, 
                                      (ctl_validation$new_batch - ctl_validation$mean)/ctl_validation$stddev, "PASS"))

# KNOWN POSITIVES  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_validation = cbind(pos, ctl_validation[rownames(pos),])
ctl_validation$cellline = sapply(strsplit(as.character(rownames(ctl_validation)), "__"), `[`, 1)
ctl_validation$antibody = sapply(strsplit(as.character(rownames(ctl_validation)), "__"), `[`, 2)

write.table(ctl_validation, file = "qc_controls.tsv", 
            sep="\t", quote = F, row.names = T, col.names = NA)

write.table(pos_validation, file = "qc_controls_positive.tsv", 
            sep="\t", quote = F, row.names = T, col.names = NA)

## QC LINEAR MODEL PLOTS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Comparing newQC with old QC with simple linear regression analysis

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
#outdat =  ctl_validation[which((ctl_validation$new_batch-ctl_validation$mean) < -twosd|(ctl_validation$new_batch-ctl_validation$mean) > twosd),]
outdat=ctl_validation
outdat$cellline = sapply(strsplit(as.character(rownames(outdat)), "__"), `[`, 1)
outdat$probe = sapply(strsplit(as.character(rownames(outdat)), "__"), `[`, 2)

baplot_ctl = ggplot(outdat) +
  geom_point(aes(x=mean, y=new_batch-mean, color=antibody)) +
  facet_wrap(~ cellline, scale="fixed") +
  geom_hline(yintercept=c(0, twosd, -twosd), color="red", linetype = 2) +
  geom_text(data=data.frame(x=0,y=c(-twosd, twosd),l=c("2SD_diff", "-2SD_diff")), aes(x=x, y=y, label=l), color=c("red"), hjust=0, size=3) +
  labs(title = "Mean-Difference Plot of Celllines", x="Mean of validation batches", y="New Batch - Mean of validation batches") + 
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size=8))

baplot_ab = ggplot(outdat) +
  geom_point(aes(x=mean, y=new_batch-mean, color=cellline)) +
  facet_wrap(~ antibody, scale="fixed") +
  geom_hline(yintercept=c(0, twosd, -twosd), color="red", linetype = 2) +
  geom_text(data=data.frame(x=0,y=c(-twosd, twosd),l=c("2SD_diff", "-2SD_diff")), aes(x=x, y=y, label=l), color=c("red"), hjust=0, size=2) +
  #geom_point(data = outdat, aes(x=mean,  y=outdat$new_batch-outdat$mean,color=outdat$probe)) +
  labs(title = "Mean-Difference Plot of Antibodies", x="Mean of validation batches", y="New Batch - Mean of validation batches") +
  theme(legend.position="right",
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size =6),
        strip.text = element_text(size=8))


print("saving QC plot")
ggsave(filename ="QC_antibody_linear_plot.pdf", device = "pdf", plot=ggReg,width = 8, height = 8)
pdf(file="QC_antibody_meandiff_plot.pdf", width=11, height=7)
print(baplot_ctl)
print(baplot_ab)
dev.off()


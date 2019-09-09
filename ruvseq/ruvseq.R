#!/usr/bin/env Rscript
usage = "
Runs Davide Risso's version of RUVseq optimized for RNAseq normalized counts
is just the weights from factor analysis removed from observed variables(genes).
Not an ideal way to perform factor analysis or perform RUV analysis.

If includeBCCL=T, reccommend count matrix, use filtperc=10, mincnt = 10, normtype=CPM, setk=7.
If includeBCCL=F, reccommend TPM matrix, use filtperc=10, mincnt = 3, normtype=NONE, setk=1.

Currently wrapped around RUVs,
Risso et.al. Nature Biotechnology volume 32, (2014)

Inputs:
1.  mat=matrix 
2.  includeBCCL = TRUE/FALSE, whether you're data includes the BCCL data 
with replicates: SkBr3, T47D, HCC1954.
3.  ctrlregex = regexpr that repeats are based off. ie. \"*UHR*\"
4.  sampregex = regexpr for samples to filter on. */*16113*/NS500642|NS500681|NB551673
5.  filtperc: minimum number of samples with counts/tpm/rpkm less that
filt_cnt filtered out.
6.  mincnt: minimum number of counts/tpm/rpkm to be filtered out
7.  normalization prior to RUV. CPM=counts per million, NONE=no normalization
8.  number of factors to remove


Outputs:
1.  normalizedCounts
2.  PCA
3.  Relative log expression plots
"

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(RUVSeq))
suppressPackageStartupMessages(library(ggplot2))

parser <- ArgumentParser()
parser$add_argument("--mat", help="table/matrix. genes are rows, samples are columns")
parser$add_argument("--includeBCCL", action="store_true", default=FALSE,
                    dest="includeBCCL", help="data includes the BCCL")
parser$add_argument("--ctrlregex", type="character", default="UHR",
                    dest="ctrlregex", help="regexpr that repeats are based off")
parser$add_argument("--sampregex", type="character", default="*",
                    dest="sampregex", help="regexpr for samples to filter on")
parser$add_argument("--filtperc", type="double", default=10,
                    dest="filtperc", help = "minimum number of samples with counts/tpm/rpkm less that filt_cnt filtered out")
parser$add_argument("--mincnt", type="double", default=3,
                    dest="mincnt", help="minimum number of counts/tpm/rpkm to be filtered out")
parser$add_argument("--normtype", type="character", default="NONE",
                    dest="normtype",help="normalization prior to RUV. CPM=counts per million, NONE=no normalization")
parser$add_argument("--setk", default=1, type="double",
                    dest="setk",help="number of factors to remove")

args <- parser$parse_args()

mat = args$mat
includeBCCL = args$includeBCCL
ctrlregex = args$ctrlregex
sampregex = args$sampregex
filtperc= args$filtperc
mincnt= args$mincnt
normtype= args$normtype
setk= args$setk

# Functions
processRUV <- function(filtmat, includeBCCL=includeBCCL, ctrlregex = ctrlregex, setk){
  #' Filter by UHR reps by default. Filtmat must include UHR replicates
  #' 
  #'  For use in normalize RUV functions
  #' @param filtmat (matrix): filtered matrix of sample=columns, genes=rows
  #' @param includeBCCL (logic): TRUE/FALSE
  #' @param setk (numeric): number of factors to remove. 
  #' if includeBCCL=T use CPM + setk=7
  #' if includeBCCL=F use TPM + setk=1
  #' @return ruvset
  #' @import RUVSeq
  #' @export
  
  filtmat= as.matrix(round(filtmat, digits=0))
  
  #set <- EDASeq::newSeqExpressionSet(as.matrix(filtmat))
  
  if (includeBCCL == TRUE){
    replicate_matx = colnames(filtmat)
    replicate_matx[grep(ctrlregex, (colnames(filtmat)), value=FALSE)] = "UHR"
    replicate_matx[grep("*SkBr3*", (colnames(filtmat)), value=FALSE)] = "SkBr3"
    replicate_matx[grep("*T47D_6*", (colnames(filtmat)), value=FALSE)]="T47D_6"
    replicate_matx[grep("*HCC1954*", (colnames(filtmat)), value=FALSE)]="HCC1954"
  } else {
    print("no BCCL")
    replicate_matx = colnames(filtmat)
    replicate_matx[grep(ctrlregex, (colnames(filtmat)), value=FALSE)] = "UHR"
  }
  scIdx = RUVSeq::makeGroups(replicate_matx)
  dim(scIdx)
  ruvs_set <- RUVSeq::RUVs(filtmat, cIdx=rownames(filtmat), k=setk, scIdx=scIdx)
  return(ruvs_set)
  
}

normalizeRUV <- function(exp_matx=exp_matx, 
                         includeBCCL=includeBCCL, 
                         sampregex=sampregex,
                         ctrlregex=ctrlregex,
                         filtperc=filtperc, 
                         mincnt=mincnt, 
                         normtype=normtype, 
                         setk) {
  #' normalize RUV
  #' 
  #' Takes a matrix and filters by counts greater than 10 in filtperc \% of SMMART samples
  #' @param exp_matx: matrix of expression, rownames=genes
  #' @param includeBCCL: TRUE OR FALSE
  #' @param sampregex: regular expression of samples to use for filtering
  #' @param ctrlregex: samples to filter the counts on
  #' @param filtperc (int): 10\% of SMMART samples will have counts greater than mincnt
  #' @param mincnt (int): TPM use 3, counts to CPM use 10
  #' @param normtype (string): "NONE", "CPM"
  #' @param setk (numeric): number of factors to remove. 
  #' if includeBCCL=T use CPM, setk=7, filtperc=10, mincnt = 10, 
  #' if includeBCCL=F use TPM, setk=1, filtperc=10, mincnt = 3, 
  #' @return list of matrices
  #' 
  #' @export
  #smmart samples differentiated by sequencer name
  samples_to_filter_on = colnames(exp_matx[,grepl(sampregex, colnames(exp_matx))])
  filtprop = as.numeric(filtperc)/100
  
  cnts = exp_matx
  
  if (normtype =="NONE"){
    #get SMMART only
    smmartcnts = cnts[,samples_to_filter_on]
    #FILTER by samples in smmart only
    perc_smmart = floor(filtprop*(ncol(smmartcnts)))
    gene_smmart = rownames(smmartcnts[apply(data.frame(smmartcnts), 1, 
                                            function(x) length(x[x>mincnt])>=perc_smmart),])
    filtmat = cnts[gene_smmart,]
    #ruv
    ruvs = processRUV(filtmat = filtmat, includeBCCL=includeBCCL, ctrlregex = ctrlregex, setk=setk)
    
  }
  
  if (normtype =="CPM"){
    #cpm 
    cpm = apply(cnts,2, function(x) (x/sum(x))*1000000)
    #get SMMART only
    smmartcpm = cpm[,samples_to_filter_on]
    #FILTER
    perc_smmart = floor(filtprop*(ncol(smmartcpm)))
    genelistcpm = rownames(smmartcpm[apply(data.frame(smmartcpm), 1, 
                                           function(x) length(x[x>mincnt])>=perc_smmart),])
    filtcpm = cpm[genelistcpm,]
    #ruv cpm
    ruvs = processRUV(filtmat = filtcpm, includeBCCL=includeBCCL, ctrlregex = ctrlregex, setk=setk)
  }
  return(ruvs)
}

# Plotting Functions

pcaplt<- function(mat, title="PCA Plot"){
  #' PCA plot function
  #' @param mat: dataframe, datatable. First column is gene names, rows are genes, columns are samples.
  #' @param title (string): Title of PCA plot.
  #' @return pcaplot
  #' @export
  lg=log(mat[,-1]+1)
  #lg=mat
  var=lg[apply(lg, 1, var, na.rm=TRUE) != 0,]
  cc.var=var[complete.cases(var),]
  
  pca_prcomp = prcomp(t(var), center = T, scale = F)
  
  PC1_and_PC2 = data.frame(PC1=pca_prcomp$x[,1], PC2= pca_prcomp$x[,2], type = rownames(pca_prcomp$x))
  PC1_and_PC2$col = ifelse(grepl(ctrlregex, PC1_and_PC2$type), "replicate", "sample")
  
  perc=(pca_prcomp$sdev^2)/sum(pca_prcomp$sdev^2) *100
  labs <- sapply(seq_along(perc), function(i) {paste("PC ", i, " (", round(perc[i], 2), "%)", sep = "")})
  
  
  p=ggplot(PC1_and_PC2,aes(PC1, PC2)) + 
    geom_point(aes(colour=col), size=4) +
    geom_text(aes(label=type, colour=col), vjust=-1) +
    scale_color_manual(values = c("sample" = "black", "replicate" = "red")) +
    labs(title=title, x=labs[1], y=labs[2]) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major=element_line(colour="gray"),
          plot.title = element_text(hjust = 0.5),
          legend.text=element_text(size=4),
          legend.position="right")
  return(p)
}


bprle <- function (mat, title="RLE", ...){
  #' Prettier relative log expression (RLE) boxplots
  #' 
  #' @param mat: make it matrix.
  #' @param title (string): title of the box plot.
  #' @return boxplot
  #' @export
  par(mar=c(20,2,1,1))
  #x=counts(set)
  y <- log2(mat + 1)
  median <- apply(y, 1, median)
  rle <- apply(y, 2, function(x) x - median)
  colors=ifelse(grepl(ctrlregex, colnames(mat)), "red", "gray")
  #boxplot(rle, ..., outline=FALSE, las=2, col=colors[pData(set)[["batch"]]], main=title, cex.axis=0.6, legend="bottom")
  boxplot(rle, ..., las=2, main=title, outline=FALSE, col=colors)
  abline(h = 0, lty = 2)
}


## MAIN

exp_matx = read.csv(mat, sep="\t", row.names = 1)



ruv_seqobj = normalizeRUV(exp_matx = exp_matx, 
                          includeBCCL = includeBCCL,
                          sampregex = sampregex,
                          ctrlregex = ctrlregex,
                          filtperc=filtperc, 
                          mincnt=mincnt, 
                          normtype=normtype, 
                          setk = setk)

norm_ruv =ruv_seqobj$normalizedCounts

pre_pca = pcaplt(exp_matx, "PCA of Input Data")
post_pca = pcaplt(norm_ruv, "PCA of RUV Processed Data")


write.table(norm_ruv, file=paste0("normRUV.txt"), sep="\t", col.names = NA, row.names = T, quote = F)
pdf("normalization_plots.pdf", pointsize=14)
pre_pca
post_pca
bprle(exp_matx, "Relative Log Exp of Input Data")
bprle(norm_ruv, "Relative Log Exp Processed Data")

dev.off()




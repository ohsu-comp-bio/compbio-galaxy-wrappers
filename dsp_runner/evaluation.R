#Evaluation Functions

#' Cluster and compute PCs on a list of data
#'
#' For each supplied matrix, the data is first centered and scaled, clusters estimated using PAM 
#' and the first two PCs recorded.
#'
#' @param k The expected number of clusters for PAM 
#' @param ... A series of numeric matrix (log2 scale) with dimensions: features (e.g. antibodies) x samples
#' with each argument named.
#' 
#' A list with elements: \describe{
#'    \item{data}{A \code{data.table} containing the sample names, the first two PCs and associated cluster assignments for each matrix}
#'    \item{avg_sil}{The average silhouette value computed for each clustering}
#' }
#' @import data.table
#' @import cluster
#' @export
run_pca <- function(k=15, ...){
  
  exprs.mats <- list(...)
  
  tmp.res <- lapply(exprs.mats, function(x){
    
    ss.x <- scale(t(x), center=T, scale=T)
    
    tmp.clusts <- pam(ss.x, k=k)
    
    si <- summary(silhouette(tmp.clusts))$avg.width
    
    tmp.p <- prcomp(ss.x, center=F, scale.=F)
    
    list(data=data.table(barcode=rownames(tmp.p$x), tmp.p$x[,1:2], clusters=tmp.clusts$cluster), avg_sil=si)
    
  })
  
  
  list(data=rbindlist(lapply(tmp.res, function(x){
    
    x$data
    
  }), idcol="type"),
  
  avg_sil=sapply(tmp.res, "[[", "avg_sil")
  
  )
  
}


.compute.tra.geo <- function(x, raw.mat, cl.repmat, abs=NULL){
  
  tmp.list <- scale_by(2^raw.mat, rows=abs, rm_cols=c(x$Var1, x$Var2))
  
  data.table(ProbeName=rownames(tmp.list$cor_mat), tra=log2(tmp.list$cor_mat)[,x$Var1] - log2(tmp.list$cor_mat)[,x$Var2], idx=x$idx, barcode=x$Var1)
  
}

.compute.tra.ruv <- function(x, raw.mat, cl.repmat, abs=NULL){
  
  tmp.repmat <- cl.repmat
  
  tmp.repmat[c(x$Var1, x$Var2),x$sample] <- 0L
  
  tmp.pcs <- compute_factors(raw.mat, tmp.repmat)
  
  tmp.norm <- apply_norm(raw.mat[,c(x$Var1, x$Var2)], tmp.pcs$pcs, k=2, controls=abs)
  
  data.table(ProbeName=rownames(tmp.norm),tra=tmp.norm[,1] - tmp.norm[,2], idx=x$idx, barcode=x$Var1)
  
}

#' Compute Technical Replicate Agreement (TRA)
#'
#' Using the replicate->sample relationships in `cl.repmat` compute the TRA defined as the pairwise
#' log fold change between replicates within a sample/type
#'
#' @param raw.mat A numeric matrix (log2 scale) with dimensions: features (e.g. antibodies) x samples
#' @param cl.repmat An indicator matrix of samples x types where a '1' indicates that a given sample Y is a replicate of type X and '0' indicates otherwise
#' @param norm.type Type of normalization, currently either Geomean or RUV
#' @param row_list A named list containing the features/antibodies to use for normalization
#' @return A \code{data.table} containing both the 'Raw' (i.e. unnormalized) and specified normalization TRA values per Probe/Ab and replicate pairs (idx).
#' @import data.table
#' 
#' @export
get_tra <- function(raw.mat, cl.repmat, norm.type=c("ruv", "geo"),row_list=list()){
  
  norm.type <- match.arg(norm.type)
  
  Var1=Var2=value=`.`=idx=NULL # due to NSE notes in R CMD check
  
  .compute.tra <- switch(norm.type, ruv=.compute.tra.ruv, geo=.compute.tra.geo)
  
  #record the possible pairs by matching on their 'parent' cell lines/tissues
  pair.dt <- data.table(reshape2::melt(cl.repmat %*% t(cl.repmat), as.is=T))[(value > 0) & (Var1 != Var2)]
  
  #as this is symmetric remove the duplicates
  pair.dt[,idx:=sapply(strsplit(paste(Var1, Var2), " "), function(x) paste(sort(x), collapse=" "))]
  pair.dt <- pair.dt[!duplicated(idx)]
  
  #add in 'tissue/cell-line' definitions
  
  cl.barcs <- data.table(reshape2::melt(cl.repmat, as.is=T))[value > 0]
  
  pair.dt <- merge(pair.dt, cl.barcs[,.(Var1, sample=Var2)], by="Var1")
  
  #for each cell line there is expected to be ncol(cl.repmat) total batches, repeated for each cell line
  #these rels won't hold if unequal number of samples per batch
  stopifnot(pair.dt[,.N] == (choose(nrow(cl.repmat)/ncol(cl.repmat), 2) * ncol(cl.repmat)))
  
  #per definition, need to re-normalize holding out each pair....
  
  split.pair <- split(pair.dt, by="idx")
  
  raw.tra <- rbindlist(lapply(split.pair, function(x){
    
    data.table(type="Raw", ProbeName=rownames(raw.mat),tra=raw.mat[,x$Var1] - raw.mat[,x$Var2], idx=x$idx, barcode=x$Var1)
    
  }))
  
  all.tra <- rbindlist(lapply(row_list, function(z){
    
    rbindlist(lapply(split.pair, .compute.tra, raw.mat=raw.mat, cl.repmat=cl.repmat, abs=z))
    
  }), idcol="type")
  
  rbind(
    raw.tra,
    all.tra
  )
  
}



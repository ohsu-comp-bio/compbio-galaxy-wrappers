#Normalization Functions

#' Perform 'Geomean' scaling normalization
#'
#'https://blog.nanostring.com/geomx-online-user-manual/Content/DataAnalysis/Algorithm_Details.htm
#'Normalization factor for each segment shall be calculated as the ratio of the average of all geometric mean or 
#'arithmetic mean of selected normalization probes across all segments to the geometric mean or arithmetic mean of 
#'selected normalization probes for that segment.
#'
#' @param mat A numeric matrix features (e.g. antibodies) x samples
#' @param rows The features to be used in the normalization, either rownames or integers
#' @param rm_cols The samples NOT used in determination of correction factor, they will still be part of the output 
#' @return A list with elements:  \describe{
#'    \item{cor_mat}{The features x samples normalized matrix}
#'    \item{cor_facs}{The computed scaling factors}
#' }
#' @export
scale_by <- function(mat, rows=NULL, rm_cols=NULL){
  
  if (missing(rows) || is.null(rows)){
    rows <- rownames(mat)
  }
  
  if (missing(rm_cols) || is.null(rm_cols)){
    keep.cols <- colnames(mat)
  }else{
    keep.cols <- setdiff(colnames(mat), rm_cols)
  }
  
  all.geo <- (exp(colMeans(log(mat[rows,,drop=F]))))
  
  all.cors <- mean(all.geo[keep.cols]) / all.geo
  
  cors.mat <- matrix(rep(all.cors, each=nrow(mat)), ncol=length(all.cors), byrow=F)
  
  res.mat <- mat * cors.mat
  
  list(cor_mat=res.mat, cor_facs=all.cors)
  
}

#' Estimate RUVIII normalization factors
#'
#' Part 1 of a re-implementation of the RUVIII function in the \code{ruv} package.
#'
#' @param mat A numeric matrix (log2 scale) with dimensions: features (e.g. antibodies) x samples
#' @param repmat An indicator matrix of samples x types where a '1' indicates that a given sample Y is a replicate of type X and '0' indicates otherwise
#' @return A list with elements: \describe{
#'    \item{fit}{A 'limma::MArrayLM' object of the model \eqn{mat ~ 0 + type1 + type2 + ... + typeN}}
#'    \item{resid}{A matrix of residuals from the model fit}
#'    \item{pcs}{A 'prcomp' object of the PCs from the model residuals}
#' }
#' @import limma
#' @import stats
#' @export
compute_factors <- function(mat, repmat){
  
  fit <- lmFit(mat[,rownames(repmat)], design=repmat)
  resid <- mat[,rownames(repmat)] - (fit$coefficients %*% t(fit$design))
  pcs <- prcomp(resid, center=F, scale.=F)
  
  list(pcs=pcs, fit=fit, resid=resid)
}

#' Apply estimated RUVIII normalization factors to a matrix
#'
#' Part 2 of a re-implementation of the RUVIII function in the \link[ruv]{RUVIII} package.
#'
#' @param mat A numeric matrix (log2 scale) with dimensions: features (e.g. antibodies) x samples
#' @param pcs An object of class \code{prcomp} generated from the \code{compute_factors} function
#' @param k The number of PC's used as part of the correction
#' @param controls The specific features used for the correction, defaults to all features
#' @return A normalized matrix of dimension features x samples
#' @import limma
#' @export
apply_norm <- function(mat, pcs, k=2, controls=NULL){
  
  if (missing(controls) || all(is.na(controls))){
    controls <- rownames(pcs$x)
  }
  
  W.fit <- lmFit(t(mat)[, controls], design=pcs$x[controls,seq_len(k)])
  new.resp <- t(mat) - W.fit$coefficients %*% t(pcs$x)[seq_len(k),]
  
  t(new.resp)
}

#' IgG background correction
#' 
#' This method implements two methods of background correction as well as the default 
#' pass-through method.  The 'log2_ratio' method simply computes the difference between
#' the experimental value and IgG assuming the data is already been log2 transformed.  The
#' 'diff' method first computes the difference between experimental values and IgG on the original
#' scale.  We then set negative values to zero and add one to every entry.  Finally, a log2 transform
#' is applied.
#' 
#' @param mat A numeric matrix (log2 scale) with dimensions: features (e.g. antibodies) x samples
#' @param igg.map A \code{data.table} containing the mapping between 'ProbeName'
#' and corresponding 'igg'.
#' @param bg.method, choice of background correction method (defaults to no correction)
#' @return A new features x samples numeric matrix after application of background correction
#' @export

bg_correct <- function(mat, igg.map=NULL, bg.method=c('none', 'log2_ratio', 'diff')){
  
  if (bg.method == 'none'){
    
    return(mat)
    
  }else if(all(c('ProbeName', 'igg') %in% names(igg.map))==F){
    
    stop("ERROR: 'igg.map' needs to have 'ProbeName' and 'igg' columns defined")
    
  }
  
  if (length(setdiff(rownames(mat), igg.map$ProbeName)) > 0){
    warning("Some ProbeNames are not defined in igg.map")
  }
  
  if (bg.method == 'log2_ratio'){
    
    igg.abund <- (mat[igg.map$ProbeName,]) - (mat[igg.map$igg,])
    
  }else if (bg.method == 'diff'){
    
    igg.abund <- (2^mat[igg.map$ProbeName,]) - (2^mat[igg.map$igg,])
    igg.abund[igg.abund < 0] <- 0
    igg.abund <- log2(igg.abund + 1)
    
  }else{
    stop("ERROR: Undefined background correction method")
  }
  
  return(igg.abund)
}


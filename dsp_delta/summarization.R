#Summarization Functions

#' Generate quantile interpolation functions from a reference set
#'
#' @param mat A numeric matrix of features (e.g. antibodies) x samples, expected to be normalized.  
#' Should be from a reference set
#' @return A named list per feature containing an \code{approxfun} function
#' @import data.table
#' @import stats
#' @export
quant_func <- function(mat){
  
  ab.interp <- lapply(setNames(rownames(mat), rownames(mat)), function(x){
    
    q.map <- data.table(v=quantile(mat[x,], seq(0, 1, by=.01)), q= seq(0, 1, by=.01))
    
    fc <- approxfun(x=q.map$v, y=q.map$q, rule=2)
    
    fc
    
  })
  
  ab.interp
}

#' Interpolate corresponding quantiles from a reference distribution
#'
#' @param interp.list A list of \code{approxfun}'s as derived from \code{quant_func}  
#' @param mat A numeric matrix of features (e.g. antibodies) x samples, expected to be normalized.  
#' @return A \code{data.table} containing the corresponding normalized data and quantile for each feature and barcode.
#' @import data.table
#' @export
get_quants <- function(interp.list, mat){
  
  Var1=Var2=value=`.`=NULL # due to NSE notes in R CMD check
  
  if (is.matrix(mat)==F){
    mat <- as.matrix(mat)
  }
  
  tmp.quants <- mapply(function(x,y){
    interp.list[[x]](y)
  }, rownames(mat), data.frame(t(mat)))
  
  if (is.matrix(tmp.quants)==F){
    ref.quants <- as.matrix(tmp.quants)
  }else{
    ref.quants <- t(tmp.quants)
  }
  
  colnames(ref.quants) <- colnames(mat)
  
  ref.melt <- data.table(reshape2::melt(mat, as.is=T))
  ref.melt <- merge(ref.melt, data.table(reshape2::melt(ref.quants, as.is=T))[,.(Var1, Var2, quant=value)], by=c("Var1", "Var2"))
  names(ref.melt) <- c("ProbeName", "barcode", "norm", "quant")
  
  ref.melt
  
}

#' Compute mean and standard deviation for a reference set
#'
#' @param mat A numeric matrix of features (e.g. antibodies) x samples, expected to be normalized.  
#' Should be from a reference set
#' @return A \code{data.table} containing the estimated mean and standard deviation from the matrix for each feature
#' @import data.table
#' @import matrixStats
#' @export
learn_mean_sd <- function(mat){
  
  data.table(ProbeName=rownames(mat), ref_mean=matrixStats::rowMeans2(mat), ref_sd=matrixStats::rowSds(mat))
  
}

#' Compute maximum correlation amonst ROIs
#'
#' @param mat A numeric matrix features (e.g. antibodies) x samples, expected to be un-normalized and log2 transformed
#' @param meta A \code{data.table} containing metadata as returned by \code{process_batches}
#' @param roi.thresh The lowest correlation value before considering the ROI(s) an outlier
#' @return A new copy of 'meta' containing the maximum correlation amongst ROIs per sample and 'rm_croi'
#' a flag indicating if the max value was less than the specified threshold
#' @import data.table
#' @import matrixStats
#' @import stats
#' @export
flag_roi <- function(mat, meta, roi.thresh=.9){
  
  tmp.meta <- copy(meta)
  
  rm_croi=max_cor=`Segment (Name/ Label)`=NULL # due to NSE notes in R CMD check
  
  melt.exp <- data.table(reshape2::melt(mat, as.is=T))
  
  melt.exp.m <- merge(tmp.meta, melt.exp, by.x="barcode", by.y="Var2")
  
  rep.cors <- rbindlist(lapply(split(melt.exp.m, by=c("Segment (Name/ Label)","sample_id", "batch")), function(x){
    
    tmp.cor <- cor(reshape2::acast(Var1~croi,value.var="value", data=x))
    
    diag(tmp.cor) <- NA
    
    data.table(`Segment (Name/ Label)`=x[1,`Segment (Name/ Label)`],sample_id=x$sample_id[1], batch=x$batch[1],croi=colnames(tmp.cor), max_cor=matrixStats::colMaxs(tmp.cor, na.rm=T))
    
  }))
  
  rep.cors[is.infinite(max_cor), max_cor:=1]
  
  if (rep.cors[,.N] != tmp.meta[,.N]){
    stop("ERROR: A size mismatch has occured, please check metadata to ensure 'Segment (Name/ Label)', 'sample_id' and 'batch' look sane")
  }
  
  tmp.meta <- merge(tmp.meta, rep.cors, by=c("Segment (Name/ Label)", "sample_id", "batch", "croi"))
  
  tmp.meta[,rm_croi:=max_cor < roi.thresh]
  
  tmp.meta
}


.summarize_roi <- function(norm.mat, meta, num.roi.avg=3){
  
  # due to NSE notes in R CMD check
  num_batch=sample_id=`Segment (Name/ Label)`=croi=use_roi=avg_barcode=`.`=NULL 
  
  ##we will first order by increasing ROI
  ord.meta <- meta[order(num_batch, sample_id, `Segment (Name/ Label)`, as.numeric(croi))]
  
  ord.meta[,use_roi:=ifelse(seq_len(.N) %in% seq_len(num.roi.avg), "avg", croi),by=.(num_batch, sample_id, `Segment (Name/ Label)`)]
  
  segment.proc <- lapply(split(ord.meta, by="Segment (Name/ Label)"), function(x){
    
    #average over the ROIs
    avg.abund <- sapply(split(x, by=c("sample_id", "num_batch", "use_roi")), function(y){
      rowMeans(norm.mat[,y$barcode,drop=F])
    })
    
    new.x <- copy(x)
    
    new.x[,avg_barcode:=paste(sample_id, num_batch, use_roi, sep=".")]
    
    list(meta=new.x[,.(num_ROI=.N),by=.(avg_barcode, sample_id=sample_id, num_batch)], avg_abund=avg.abund)
    
  })
  
  segment.proc
  
}

#' High-level procedure for pre-processing DSP data using TMAs
#'
#' This function first computes RUVIII normalization factors after background correction
#' for the tumor microarrays (TMAs).  These factors are then applied to the experimental
#' data after background correction.  The resulting matrix is then averaged over the ROIs
#' to produce normalized antibody x sample matrices for each segment.
#'
#' @param tma.meta A \code{data.table} containing the following columns: \describe{
#'    \item{num_batch}{Numeric batch identifier or other ordered run signifier such as date}
#'    \item{barcode}{The unique sample identifier for the TMA sample/run}
#'    \item{name}{Harmonized TMA sample name, consistent across batches}
#'    \item{type}{The 'Type' label of TMA sample, only sample's with values in 
#'    \code{use.type}will be used}
#' }
#' @param tma.mat A numeric matrix (log2 scale) with dimensions: features (e.g. antibodies) x TMA sample barcodes
#' @param exp.meta A \code{data.table} needs to have the following columns: \describe{
#'    \item{Segment (Name/ Label)}{Segment Label}
#'    \item{sample_id}{Sample Identifier}
#'    \item{barcode}{The unique sample identifier for the sample/run}
#'    \item{num_batch}{Numeric batch identifier or other ordered run signifier such as date}
#'    \item{croi}{Corrected region of interest typically as generated in a previous step}
#' }
#' @param exp.mat A numeric matrix (log2 scale) containing experimental data with 
#' dimensions: features (e.g. antibodies) x (# segments x # rois x # samples)
#' @param igg.map A \code{data.table} containing the mapping between 'ProbeName'
#' and corresponding 'igg'.
#' @param bg.method, choice of background correction method (defaults to no correction)
#' @param controls The specific features used for the correction, defaults to all features
#' @param use.type A character vector of values in the \code{tma.meta} \code{type} column to be used.
#' Defaults to 'quant'.
#' @param k The number of PC's used as part of the correction
#' @param num.roi.avg The number of ROIs to average, assumes they are numbered as the first 1:num.roi.avg
#' @return A list with one element per segment each of which contains: \describe{
#'    \item{meta}{A summarized meta \code{data.table} containing the number of ROIs, the new sample barcodes ('avg_barcode'), sample ID and averaged abundance }
#'    \item{avg_abund}{A features (e.g. antibodies) x samples numeric matrix}}
#'    
#' @import data.table
#' @export
preprocess_dsp_tma <- function(tma.meta, tma.mat, exp.meta, exp.mat, igg.map=NULL, bg.method=c('none', 'log2_ratio', 'diff'), controls=NULL, use.type='quant', k=2, num.roi.avg=3){
  
  tmpval=type=sample_id=num_batch=NULL # due to NSE notes in R CMD check
  
  #compute normalization factors from tma
  
  meta.cp <- copy(tma.meta)
  meta.cp[,tmpval:=1]
  cl.repmat <- reshape2::acast(barcode~name, value.var="tmpval", data=meta.cp[type %in% use.type], fun.aggregate = function(x) as.integer(length(x) > 0))
  
  if (all(colSums(cl.repmat) != meta.cp[,length(unique(num_batch))])){
    stop("ERROR: Samples are not found for every batch")
  }
  
  bg.tma.mat <- bg_correct(tma.mat, igg.map, bg.method)
  
  cl.pcs <- compute_factors(bg.tma.mat, cl.repmat)
  
  #apply bg correction to experimental and apply normalization
  
  bg.exp.mat <- bg_correct(exp.mat, igg.map, bg.method)
  
  #apply RUV normalization, note that it is technically applied sample-by-sample so ok to mix tumor/stroma segments here
  norm.igg.abund  <- apply_norm(bg.exp.mat, cl.pcs$pcs, k, controls)
  
  .summarize_roi(norm.igg.abund, exp.meta, num.roi.avg)
  
}

#' Cohort-level scoring of abundance data
#' 
#' Scores antibofy abundance data with respect to a reference cohort returning either the estimated reference quantile ('quant') or
#' a robust version of the Zscore ('rscore') indicating deviation from the reference cohort median.  If the set of reference samples
#' isn't supplied, uses the entire set of input samples as the reference cohort therefore producing an intra-cohort scoring.
#' 
#' @param norm.list A list per segment as generated by `preprocess_dsp_tma` with each element 
#' containing a list with two elements: \describe{
#'    \item{meta}{A summarized meta \code{data.table} containing the number of ROIs, the new sample barcodes ('avg_barcode'), sample ID and averaged abundance }
#'    \item{avg_abund}{A features (e.g. antibodies) x samples numeric matrix}
#'    }
#' @param ref.samples A charancter vector of barcodes corresponding to `avg_barcode` containing the samples to use as the reference cohort. 
#' If missing or NULL will treat every sample as part of the reference. 
#' @param score.type One of either reference quantile ('quant') or robust Zscore ('rscore')
#' @return A list with one element per segment: \describe{
#'    \item{scores}{A \code{data.table} containing the antibody name (ProbeName), barcode, normalized abundance and a corresponding 'quant' or 'rscore' column }
#'    \item{ref_abund}{A \code{data.table} containing the antibody name (ProbeName), barcode and normalized abundance for the reference cohort}
#'    }
#' @import matrixStats
#' @import data.table
#' @export
score_abs <- function(norm.list, ref.samples=NULL, stroma=T, score.type=c("quant", "rscore")){

  score.type <- match.arg(score.type)

  mads=sample_id=avg_barcode=rscore=medians=score_abs=NULL # due to NSE notes in R CMD check

  no.ref <- F

  if (missing(ref.samples) || is.null(ref.samples)){
    no.ref <- T
  }

  if (stroma){
    lapply(norm.list, function(nl){

      if (no.ref){
        tmp.ref <- nl$avg_abund
        tmp.exp <- nl$avg_abund
      }else{
        tmp.ref <- nl$avg_abund[,nl$meta[sample_id %in% ref.samples,avg_barcode],drop=F]
        tmp.exp <- nl$avg_abund[,nl$meta[sample_id %in% ref.samples == F,avg_barcode],drop=F]
      }

      if (score.type == "quant"){

        tmp.interp <- quant_func(tmp.ref)

        tmp.quants <- get_quants(tmp.interp, tmp.exp)
        names(tmp.quants)[2] <- "avg_barcode"

      }else{

        #add in robust zscore as an alternative
        #https://stats.stackexchange.com/questions/523865/calculating-robust-z-scores-with-median-and-mad

        ref.mstats <- data.table(ProbeName=rownames(tmp.ref), mads=rowMads(tmp.ref, constant=1), medians=rowMedians(tmp.ref))

        tmp.quants <- setNames(data.table(reshape2::melt(tmp.exp, as.is=T)), c("ProbeName", "avg_barcode","norm"))
        tmp.quants <- merge(tmp.quants, ref.mstats, by="ProbeName")
        tmp.quants[,rscore:=(norm-medians)/mads]
        tmp.quants[,`:=`( mads=NULL, medians=NULL)]

      }

      list(
        scores=tmp.quants,
        ref_abund=setNames(data.table(reshape2::melt(tmp.ref, as.is=T)), c("ProbeName", "avg_barcode","norm"))
      )

    })
  }else{
    nl <- norm.list
    if (no.ref){
      tmp.ref <- nl$avg_abund
      tmp.exp <- nl$avg_abund
    }else{
      tmp.ref <- nl$avg_abund[,nl$meta[sample_id %in% ref.samples,avg_barcode],drop=F]
      tmp.exp <- nl$avg_abund[,nl$meta[sample_id %in% ref.samples == F,avg_barcode],drop=F]
    }

    if (score.type == "quant"){

      tmp.interp <- quant_func(tmp.ref)

      tmp.quants <- get_quants(tmp.interp, tmp.exp)
      names(tmp.quants)[2] <- "avg_barcode"

    }else{

      #add in robust zscore as an alternative
      #https://stats.stackexchange.com/questions/523865/calculating-robust-z-scores-with-median-and-mad

      ref.mstats <- data.table(ProbeName=rownames(tmp.ref), mads=rowMads(tmp.ref, constant=1), medians=rowMedians(tmp.ref))

      tmp.quants <- setNames(data.table(reshape2::melt(tmp.exp, as.is=T)), c("ProbeName", "avg_barcode","norm"))
      tmp.quants <- merge(tmp.quants, ref.mstats, by="ProbeName")
      tmp.quants[,rscore:=(norm-medians)/mads]
      tmp.quants[,`:=`( mads=NULL, medians=NULL)]

    }

    list(
      scores=tmp.quants,
      ref_abund=setNames(data.table(reshape2::melt(tmp.ref, as.is=T)), c("ProbeName", "avg_barcode","norm"))
    )

  }
}

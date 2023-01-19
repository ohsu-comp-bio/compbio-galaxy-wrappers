#' Produce a Pathway Summary Plot
#'
#' For each patient sample plot the reference quantile per antibody along with the median per pathway (pathway score).  If multiple samples exist, then also
#' compute the difference between the last and first patient timepoint as defined as well as produce the difference between the pathway scores (difference of medians).
#' If data for a reference cohort is also provided, will draw also show boxplots of reference distribution.
#'
#' @param score.dt A \code{data.table} containing at least the following columns:\describe{
#'    \item{avg_barcode}{The unique sample barcode for the ROI averaged data}
#'    \item{sample_ord}{An ordered factor containing the ordered relationships of the patient samples (e.g. pre-dx, treatment 1, treatment 2)}
#'    \item{segment_label}{The human interpretable labels for the segments}
#'    \item{patient_ord}{An ordered factor containing the desired patient labels }
#'    \item{path_ord}{An ordered factor containing the pathway order}
#'    \item{ab_ord}{An ordered factor containing the antibody order}
#'    \item{quant}{The observed quantile relative to a reference distribution}
#'    \item{norm}{The normalized abundance values}}
#' @param ref.dt An optional \code{data.table} containing at least the following columns:\describe{
#'    \item{segment_label}{The human interpretable labels for the segments}
#'    \item{path_ord}{An ordered factor containing the pathway order}
#'    \item{ab_ord}{An ordered factor containing the antibody order}
#'    \item{norm}{The normalized abundance values}}
#' @param no.median A character vector containing the names of the \code{path_ord} levels to not summarize into pathway scores but instead leave
#' at the antibody level.
#' 
#' 
#' @return A list of lists containing ggplot2 or patchwork'd objects, with the outer list indexing segment and the inner list indexing patient
#' @import data.table
#' @import ggplot2
#' @import patchwork
#' @import stats
#' @import scales
#' @export
loli_plot <- function(score.dt, ref.dt=NULL, no.median=c("Tumor Markers", "Other Markers", "IO CLIA", "IO RUO")){

  avg_barcode=`.`=path_ord=quant=sample_ord=ab_ord=med=segment_label=patient_ord=N=med_diff=NULL # due to NSE notes in R CMD check
  
  tmp.score <- copy(score.dt) 
  
  should.box <- F
  
  if (missing(ref.dt)==F && is.null(ref.dt)==F){
    tmp.ref <- copy(ref.dt)
    should.box <- T
  }
  
  lapply(split(tmp.score, by="segment_label"), function(seg.dt){
    
    lapply(split(seg.dt, by="patient_ord"), function(samp.dt){
      
      meds <- samp.dt[as.character(path_ord) %in% no.median==F,.(med=median(quant)),by=.(path_ord, sample_ord)]
      
      unique.barcs <- samp.dt[,.N,by=avg_barcode]
      unique.barcs[,col:=scales::hue_pal()(.N)]
      
       if (should.box){
        box.plot <- ggplot(data=tmp.ref[segment_label == samp.dt$segment_label[1]], 
                           mapping=aes(y=ab_ord, x=norm)) + 
          geom_boxplot(outlier.shape=NA) + 
          geom_point(shape=21) +
          geom_vline(xintercept=0, linetype="dashed") +
          geom_point(data=samp.dt, mapping=aes(fill=avg_barcode), size=2, shape=21) +
          scale_fill_manual(values=setNames(unique.barcs$col, unique.barcs$avg_barcode), guide="none") +
          facet_grid(path_ord~segment_label, scales="free", space="free_y", labeller=.uc.labs) +
          theme_bw() + xlab("Abundance") + ylab("") +  
#          theme(axis.text.y=element_text(size=14), strip.background.y = element_blank(), strip.text.y = element_blank()) +
          ggtitle(samp.dt[1,paste(patient_ord, '|', num_batch)])

         }

      tmp.plot <- ggplot(data=samp.dt, mapping=aes(y=ab_ord, x=quant, fill=avg_barcode)) + 
        geom_segment(mapping=aes(x=0, y=ab_ord, xend=quant, yend=ab_ord), size=1.5) +
        geom_point(size=4, shape=21) + 
        scale_fill_manual(values=setNames(unique.barcs$col, unique.barcs$avg_barcode), guide="none") +
        geom_vline(data=meds, mapping=aes(xintercept=med), linetype="dashed") +
        facet_grid(path_ord~sample_ord, scales="free", space="free_y", labeller=.uc.labs) + 
        theme_bw() + ylab("") + xlab("Reference Quantile") +
        ggtitle(samp.dt[1,paste(.uc.labs(segment_label), "Segment", patient_ord)])
      
      if (should.box){
        
        tmp.plot <- tmp.plot + theme(
          axis.ticks.y = element_blank(),
          axis.text.y=element_blank()
        )
        
      }else{
        
        tmp.plot <- tmp.plot + theme(
          axis.text.y=element_text(size=14)
        )
        
      }
      
      #determine if we can add a 'delta' comparing the last value to the first
      samp.comps <- unique(samp.dt[,.(sample_ord)])[order(sample_ord)]
      
      if (samp.comps[,.N] > 1){
        
        tmp.plot <- tmp.plot +
          theme(
            strip.background.y = element_blank(),
            strip.text.y = element_blank()
          )
        
        #delta is computed last - first
        samp.delta <- samp.dt[sample_ord %in% samp.comps[c(1, nrow(samp.comps))]$sample_ord, .(med_diff=.SD[order(sample_ord),diff(quant)], .N) , by=.(ab_ord, path_ord)]
        
        stopifnot(samp.delta[,all(N==2)])
        
        delta.str <-  paste(.uc.labs(rev(samp.comps[c(1, nrow(samp.comps))]$sample_ord)), collapse=" - ")
        
        meds.delta <- meds[order(sample_ord),.(med=diff(med)), by=.(path_ord)][,.(path_ord, sample_ord=paste("Delta (",delta.str,")"), med)]
        
        path.plot <- ggplot(data=samp.delta, mapping=aes(y=ab_ord, x=med_diff, fill=med_diff)) + 
          geom_segment(mapping=aes(x=0, y=ab_ord, xend=med_diff, yend=ab_ord), size=1.5) +
          geom_point(size=4, shape=21) + 
          scale_fill_gradient2(low="blue", mid="white", high="red", guide="none") +
          geom_vline(data=meds.delta, mapping=aes(xintercept=med), linetype="dashed") +
          #scale_color_manual(values=c(increasing="red", neither="black", decreasing="blue"), guide="none") +
          scale_x_continuous(breaks=c(-1, -.5, -.1, 0, .1, .5, 1)) +
          facet_grid(path_ord~sample_ord, scales="free", space="free_y", labeller=.uc.labs) + theme_bw() +
          theme(axis.ticks.y = element_blank(),
                axis.text.y=element_blank()) + xlab("") + ylab("")
        
        if(should.box){
          
          box.plot + tmp.plot + path.plot + plot_layout(ncol=3, widths=c(1, 1.5,1))
          
        }else{
          
          tmp.plot + path.plot + plot_layout(widths=c(1.5,1))
        }
        
      }else{
        
        if(should.box){

          box.plot <- box.plot + theme(axis.text.y=element_text(size=14))
          box.plot + plot_layout(ncol=1, widths=c(1, 1.5))
          
        }else{
          tmp.plot
        }
      }
      
    })
    
  })
  
}


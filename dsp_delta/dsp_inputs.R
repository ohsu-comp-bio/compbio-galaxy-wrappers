.parse.table <- function(dt){
  
  Segment.display.name=NULL# due to NSE notes in R CMD check
  
  # JHL: Changed to Segment.display.name
  hash.locs <- dt[,grep("#",  Segment.display.name)]
  
  meta <- rbind(
    dt[1:(hash.locs[1]-1)],
    dt[(hash.locs[1]+1):(hash.locs[2]-1)]
  )
  
  p.meta <- t(meta)
  p.meta <- p.meta[grepl("^X\\d+",rownames(p.meta))==F,]
  
  meta.head <- p.meta[1,]
  p.meta <- p.meta[-1,]
  
  p.meta.dt <- data.table(barcode=rownames(p.meta),p.meta)
  names(p.meta.dt)[2:ncol(p.meta.dt)] <- meta.head
  
  count.dt <- dt[(hash.locs[2]+1):dt[,.N]]
  names(count.dt)[1:4] <- sub("\\s+", "", capwords(sub("#", "", sub("\\s+\\(.+", "", unname(unlist(dt[hash.locs[2],1:4,with=F]))))))
  
  #new version lists TargetName instead of ProbeName
  names(count.dt)[names(count.dt) == "TargetName"] <- "ProbeName"
  names(count.dt)[names(count.dt) == "TargetGroup"] <- "ProbeGroup"
  
  #also coerce numeric-like cols to numeric
  
  list(abundance=count.dt[,lapply(.SD, .check.numeric)], meta=p.meta.dt[,lapply(.SD, .check.numeric)])
}

#' Read in and pre-process Excel output from GeoMx
#'
#' @param batch.dt A \code{data.table} with two columns \code{files} indicating the 
#' input files and \code{batch} which indicates the grouping of samples (usually corresponding to a directory).
#' There should at least be one Excel file indicated by 'Initial' in the filename which has the initial (non-QC scaled) data.
#' @param exp.regex The regular expression used to identify the files containing the 'initial' (non-qc adjusted) exported data.
#' @return A list with elements: A list with elements: \describe{
#'    \item{meta}{A \code{data.table} containing the metadata including three additional derived columns: GeoMX.ID which is derived from the
#'    first portion of barcode as well as 'croi' (corrected region of interest) and 'secondary' which are derived from the second part of the barcode.  
#'    The 'secondary' ID and the 'GeoMx.ID' should correspond to the same value.}
#'    \item{data}{A \code{data.table} containing feature metadata in rows and sample barcodes as columns with corresponding numeric values}
#' }
#' @import data.table
#' @import openxlsx
#' @export
process_batches <- function(batch.dt, exp.regex="Initial", sheet="Exported dataset"){
  
  ProbeName=tmp_roi=sample_id=barcode=NULL # due to NSE notes in R CMD check
  
  split.batch <- split(batch.dt, by="batch")
  
  res.list <- lapply(split.batch, function(x){
    
    init.file <- grepl(exp.regex,x$files)
    qc.file <- grepl("qc",x$files)
    
    stopifnot(sum(init.file)==1)

    # JHL: added parameter for sheet name
    init.res <- .parse.table(data.table(openxlsx::read.xlsx(x$files[init.file], sheet=sheet)))
    
    #initial geomean scaling relative to the HYB-POS per the QC data output by the machine
    
    tmp.qc <- cbind(init.res$abundance[,1:4, with=T],
                    scale_by(as.matrix(init.res$abundance[,5:ncol(init.res$abundance)]) ,rows=which(init.res$abundance$ProbeName == "HYB-POS"))$cor_mat
    )
    
    tmp.qc <- tmp.qc[ProbeName %in% c("HYB-POS", "HYB-NEG")==F]
    
    list(meta=init.res$meta,initial=init.res$abundance, qc=tmp.qc)
    
  })
  
  #now process the metadata
  
  #JHL: Some of these spreadsheets don't have the same "LOT" rows headers
  qc.meta <- rbindlist(lapply(res.list, function(x) x$meta), idcol="batch", fill=TRUE)
  
  #unfortunately, the safest place to get these IDs is from the barcode
  qc.meta[,tmp_roi:=sapply(strsplit(qc.meta$barcode, "\\.+\\|\\.+"), "[",2)]
  # JHL: Since sample id's are now sitting in Scan_ID, will create a tmp_roi_scan to grab sample_ids.  
  qc.meta[,tmp_roi_scan:=sapply(strsplit(qc.meta$Scan_ID, "[_ ]"), "[",1)]
  
  # JHL: Changed regular expression (underscore instead of hyphen) since these are now of form 001_Tonsil  
  qc.meta[,c("croi", "sample_id"):=data.table(t(sapply(strsplit(qc.meta$tmp_roi, "[\\.\\_]+"), function(x){
    if (length(x) == 1){
      c(x, NA_character_)
    }else{
      tmp.x <- x[order(grepl("^0\\d{2}$", x), decreasing=T)]
      c(tmp.x[1], paste(tmp.x[-1], collapse="."))
    }
    
  })))]

  # JHL: Now grabbing sample_id from tmp_roi_scan if it is NA.
  #  qc.meta[is.na(sample_id)==T,sample_id:=sapply(strsplit(sapply(strsplit(barcode, "\\.+\\|\\.+"), "[",1), "\\."), function(x) paste(x[-length(x)], collapse="."))]
  qc.meta[is.na(sample_id)==T,sample_id:=tmp_roi_scan]
  
  qc.meta[,`:=`(tmp_roi=NULL)]
  # JHL: Need to remove tmp_roi_scan also.  
  qc.meta[,`:=`(tmp_roi_scan=NULL)]
  
  list(meta=qc.meta, data=res.list)
  
}


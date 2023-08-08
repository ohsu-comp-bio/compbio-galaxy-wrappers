.parse.table <- function(dt){

  hash.locs <- dt[,grep("#",  Segment.displayed.name)]

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
#' @return A list with elements: A list with elements: \describe{
#'    \item{meta}{A \code{data.table} containing the metadata including three additional derived columns: GeoMX.ID which is derived from the
#'    first portion of barcode as well as 'croi' (corrected region of interest) and 'secondary' which are derived from the second part of the barcode.
#'    The 'secondary' ID and the 'GeoMx.ID' should correspond to the same value.}
#'    \item{data}{A \code{data.table} containing feature metadata in rows and sample barcodes as columns with corresponding numeric values}
#' }
process_batches <- function(batch.dt, exp.regex="Initial", geo_tma_lbl=TRUE){

  split.batch <- split(batch.dt, by="batch")

  res.list <- lapply(split.batch, function(x){

    init.file <- grepl(exp.regex,x$files)
    qc.file <- grepl("qc",x$files)

    stopifnot(sum(init.file)==1)

    init.res <- .parse.table(data.table(openxlsx::read.xlsx(x$files[init.file], sheet="Export dataset")))

    which.pos <- init.res$abundance[ProbeName=="HYB-POS",5:ncol(init.res$abundance)]

    #need to update this part with scale_by
    cor.fac <- unlist(which.pos) / mean(unlist(which.pos))

    tmp.cor <- matrix(rep(cor.fac, each=init.res$abundance[,.N]), ncol=length(cor.fac), byrow=F)

    tmp.qc <- cbind(init.res$abundance[,1:4, with=T], init.res$abundance[,5:ncol(init.res$abundance)] / tmp.cor)

    tmp.qc <- tmp.qc[ProbeName %in% c("HYB-POS", "HYB-NEG")==F]

    #for some, qc exports exist, check the sanity of the above code
    if (x[,.N] == 2 & sum(qc.file)==1){
      qc.res <- .parse.table(data.table(openxlsx::read.xlsx(x$files[qc.file], sheet="Export dataset")))

      stopifnot(isTRUE(all.equal(qc.res$abundance, tmp.qc)))

    }

    list(meta=init.res$meta,initial=init.res$abundance, qc=tmp.qc)

  })

  #now process the metadata

  qc.meta <- rbindlist(lapply(res.list, function(x) x$meta), idcol="batch")

  qc.meta[,tmp_roi:=sapply(strsplit(qc.meta$barcode, "\\.+\\|\\.+"), "[",2)]
  qc.meta[,c("secondary", "croi"):=data.table(t(sapply(strsplit(tmp_roi, "\\.*\\-\\.*"), function(x){
    if (length(x) == 1){
      c(NA, x)
    }else{
      x[order(grepl("^0\\d{2}$", x))]
    }

  })))]


  qc.meta[,GeoMx.ID:=sapply(strsplit(sapply(strsplit(qc.meta$barcode, "\\.+\\|\\.+"), "[",1), "\\."), function(x) paste(x[-length(x)], collapse="."))]

  qc.meta[is.na(secondary)==F,GeoMx.ID:=ifelse(str_detect(GeoMx.ID, secondary), secondary, GeoMx.ID)]

  if(geo_tma_lbl == T) {
    qc.meta[grepl("Geometric.Segment", barcode), GeoMx.ID:="TMA"]
  }

  qc.meta[,tmp_roi:=NULL]

  list(meta=qc.meta, data=res.list)

}


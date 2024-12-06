input_qc <- function(){
  
  # Check for duplicate run IDs within batch files
  if(length(unique(batch.dt$files)) != length(unique(batch.dt$batch))){
    message("Check for redundant run ID(s).")
  }
  
  # Check for duplicate TMA/samples
  dup_idx <- duplicated(qc.meta$barcode)
  dupes <- qc.meta[dup_idx]
  if(nrow(dupes)>0){
    message('Error -- Duplicate TMA/samples: ', as.vector(dupes$barcode))
    quit(status=5)
  }
  
  # Check all batches to ensure they contain at least one clinical sample
  for (i in unique(qc.meta$batch)){
    clin_only <- qc.meta %>% select(batch, ROI_ID) %>% filter(batch==i) %>% drop_na()
    if(nrow(clin_only)==0){
      message(paste0('Error -- batch ', i, ' does not contain any clinical samples.'))
    }
  }
}

remove_dummy <- function(){
  # Search for research batches
  research_list <- batch.dt[str_detect(batch.dt$files, 'Research'),]
  research_list <- as.vector(research_list$batch)
  # Remove from meta
  qc.meta <- qc.meta %>% filter(is.na(`ROI_ID`) & batch %in% research_list | !batch %in% research_list)
  qc.meta
}
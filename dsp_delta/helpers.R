#from toupper examples--not part of unit tests
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

.check.numeric <- function(x){
  
  num.x <- suppressWarnings(as.numeric(x))
  
  if (all(is.na(num.x))==F){
    num.x
  }else{
    x
  }
  
}

.uc.labs <- function(label, multi_line = TRUE){
  
  labels <- lapply(label, function(x) capwords(as.character(x)))
  labels
}


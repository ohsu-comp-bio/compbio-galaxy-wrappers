## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
listofpackages=c("knitr", "reshape2", "ggplot2", "rmarkdown")
for (i in listofpackages){
  if(! i %in% installed.packages()){
    install.packages(i, dependencies=TRUE, repos="https://cran.rstudio.com/")
  }}
lapply(listofpackages, require, character.only = TRUE)
# ---- md ----
usage = paste0("Rscript drugprediction_report.R basemat_filename response sample prediction_file drug_table_file")

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 6) {
  stop(print(usage), call.=FALSE)
} else {
  tooldir = args[1]
  basemat_filename = args[2]
  response = args[3]
  sample = args[4]
  prediction_file = args[5]
  drug_table_file = args[6]
  print(args)
  
  rmarkdown::render(paste0(tooldir,"/drugprediction_report.Rmd"),
                    params = list(basemat_filename = basemat_filename,
                                  response = response,
                                  sample = sample,
                                  prediction_file=prediction_file,
                                  drug_table_file = drug_table_file))
}

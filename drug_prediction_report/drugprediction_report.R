#!/usr/bin/env R

library(knitr)
library(reshape2)
library(kableExtra)
library(ggplot2)



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

    rmarkdown::render(paste0(tooldir,"/drugprediction_report.Rmd"), params = list(basemat_filename = basemat_filename,
                                                 response = response,
                                                 sample = sample,
                                                 prediction_file=prediction_file,
                                                 drug_table_file = drug_table_file))

}


#!/usr/bin/env Rscript

## Finds output files for output directory and makes report
## 
##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(pdftools))
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

version="0.1"
mat_version = "20200512"
ihc_version = "20200507"

parser <- ArgumentParser()
parser$add_argument("--repo", type="character", default="/Users/patterja/Workspace/nanostring/nanostring", 
                    dest="repo", help="code directory name. Needed to find Rmd file.")
parser$add_argument("--rawdata", help="rawdata.txt", dest="input_file")
parser$add_argument("--md_file", type="character", default= "/Users/patterja/Box/NANOSTRING/nanostring_metadata.xlsx",
                    dest="mbc_md_file", help="metastatic breast cancer metadata file")
parser$add_argument("--version", action="version", version=paste0('%(prog)s = ', version))

args <- parser$parse_args()
repo = args$repo
rawdata = args$input_file
md_file = args$mbc_md_file


outdir = dirname(normalizePath(rawdata))
batch=basename(outdir)

#Things to include and exclue
controls=c("MCF7","HCC1954","BT474","HeyA8","MDA468.control","MDA468+EGF")
omit = c("Histone H3", "S6", "RbAb-IgG", "MmAb-IgG1", "p-TSC2", "TSC2")
omitregex = paste0(paste0("^", omit), collapse = "|")
ab.ctrl = "IgG|POS|NEG|^S6|^Histone"

# metadata
md = read.xlsx(xlsxFile=md_file, sheet="nanostring_metadata", check.names=T)
md$sampcolumn = make.names(paste0(md$Batch, "__", md$Sample.Name))

# NEW BATCH
new_batch = read.table(file = rawdata, sep="\t", row.names=2, stringsAsFactors=F, header=T, check.names = T)
new_batch[,c("CodeClass", "Accession")] <- NULL

# loop over samples
for (samp in setdiff(colnames(new_batch), make.names(controls))){
  print(samp)
  if (grepl("blank",samp,ignore.case=TRUE)){
    print("sample blank")
  } else{
    rmarkdown::render(file.path(repo, "nanostring_report.Rmd"), 
                      output_dir = outdir,
                      output_file = sprintf("%s_nanostring_report.pdf",samp),
                      params = list(samp = samp, batch = batch, data_dir=outdir))
    
  }
}


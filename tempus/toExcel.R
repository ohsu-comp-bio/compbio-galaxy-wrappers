

xlsxExists <- require("xlsx")

if (!xlsxExists){
	install.packages("xlsx")
	library("xlsx")
}

argparseExists <- require("argparse")

if (!argparseExists){
	install.packages("argparse")
	library("argparse")
}

parser <- ArgumentParser()

parser$add_argument("--tsv", nargs="*")
parser$add_argument("--png", nargs="*")
args <- parser$parse_args()


workbook <- createWorkbook(type="xlsx")
COUNT <- 0
for (file in args$png){
	COUNT <- COUNT + 1
	sheet <-createSheet(workbook, sheetName=paste0("Plot", COUNT))
	addPicture(file, sheet, scale = 1, startRow = 4, startColumn = 1)
}

saveWorkbook(workbook, "excelOutput.xlsx")

COUNT <- 0
for (file in args$tsv){
	COUNT <- COUNT + 1
	df <- read.table(file, header = TRUE, sep="\t")
	write.xlsx(x=df, file="excelOutput.xlsx", sheetName=paste0("Table", COUNT), append=TRUE)
}


load("mnp.RData")

library(devtools)

install_github("mwsill/mgmtstp27",subdir="/trunk/Rpackage/mgmtstp27")

install.packages("mnp.v11b4_0.1.124.tar.gz", repos = NULL, type="source", verbose=TRUE)

install.packages("mnpqc_0.1.0.tar.gz", repos = NULL, type="source", verbose=TRUE)

library(mnp.v11b4)
library(mnpqc)

args <- commandArgs(TRUE)

package.version("mnpqc")
package.version("mnp.v11b4")
path <- args[1]

sample <- args[2]

sample_path <- paste(path,sample,sep='/')

# read data using standard minfi method
RGset <- read.metharray(sample_path, verbose=FALSE)
RGset

system.time(
MNPreport(RGset,sampleID=sample,tsne = TRUE)
)
sessionInfo()

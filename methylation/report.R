
load("mnp.RData")

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

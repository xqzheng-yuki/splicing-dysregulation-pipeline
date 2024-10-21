
library(dplyr, lib.loc = "/usr/local/lib/R/site-library")
library(gplots, lib.loc = "/usr/local/lib/R/site-library")
library("rtracklayer")

# access datadetach("package:here", unload = TRUE)
setwd("/mnt/gtklab01/xiaoqing/star")

# fetch GTF
GRCh <- import("/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf")
# Filter for standard chromosomes (?)
# modify chromosome names (?)
## check naming style for gtf, bam, and bw
seqlevelsStyle(GRCh)
bamfile <- import("./unmappedAligned.sortedByCoord.out_CTX_104.bam")
bw <- import("./unmapped_CTX_104.bw")
seqlevels(bw)
seqlevels(GRCh)
seqlevels(bamfile)

# Get bigwig
CTX104 <- BigWigFile("./unmapped_CTX_104.bw")
CTX104GR <- import(CTX104)

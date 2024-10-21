library(dplyr, lib.loc = "/usr/local/lib/R/site-library")
library(gplots, lib.loc = "/usr/local/lib/R/site-library")
library(rtracklayer)
library(GenomicFeatures)
library(Gviz)


setwd("/mnt/gtklab01/xiaoqing/star")

# fetch GTF
GRCh <- import("/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf")
GRCh_S <- keepStandardChromosomes(GRCh,pruning.mode="coarse")
txdb <- makeTxDbFromGRanges(GRCh)
txdb_S <- makeTxDbFromGRanges(GRCh_S)
keytypes(txdb)
columns(txdb)
# modify chromosome names (?)
## test-check naming style for gtf, bam, and bw
# seqlevelsStyle(GRCh)
# bamfile <- import("./unmappedAligned.sortedByCoord.out_CTX_104.bam")
# bw <- import("./unmapped_CTX_104.bw")
# seqlevels(bw)
# seqlevels(GRCh)
# seqlevels(bamfile)
# Filter for standard chromosomes (?)

### Get bigwig
import_bigwig <- function(file) {
  bw_file <- BigWigFile(file)  # Create a BigWigFile object
  data<-import(bw_file)  # Import the data into a GRanges object
  keepStandardChromosomes(data,pruning.mode = "coarse")
}
# get all bigwig file in path
bw_files <- list.files(path = ".", pattern = "^unmapped_CTX_([0-9]+)\\.bw$", full.names = TRUE)
# grouping import
tags <- as.numeric(gsub("^unmapped_CTX_([0-9]+)\\.bw$", "\\1", basename(bw_files)))
treatment_group <- c(104,108,128,154)
control_group <- c(120,125,147,148)
treatment_files <- bw_files[tags %in% treatment_group]
control_files <- bw_files[tags %in% control_group]
treatment_list <- lapply(treatment_files, import_bigwig)
control_list <- lapply(control_files, import_bigwig)
names(treatment_list) <- basename(treatment_files)
names(control_list) <- basename(control_files)
# separate import
for (i in seq_along(bw_files)){
  tag <- tags[i]
  var_name <- paste0("CTX_",tag)
  assign(var_name, import_bigwig(bw_files[i]))
}
CTX104_C <- keepStandardChromosomes(CTX_104,pruning.mode = "coarse")
seqlevels(CTX104_C)
# goi=Adnp2, (chr12:56,694,976-56,714,605)
chr_no <- "chr18"
chr_start <- 80126311
chr_end <- 80151482
goi <- GRanges(seqnames = chr_no, ranges = IRanges(start=chr_start,end=chr_end),strand='-')
subsetByOverlaps(CTX_104,goi)
# VIS
g_axis <- GenomeAxisTrack()
# make transcripts track
grTrack <- GeneRegionTrack(txdb_S, 
                           chromosome = chr_no, # chromosome number
                           start=chr_start, # start of region 
                           end=chr_end, 
                           cex.group=0.4, # display the labels a bit smaller 
                           showId=TRUE)
plotTracks(c(g_axis,grTrack))
# make bigwig data track
TreatmentTrack <- DataTrack(range = treatment_list, 
                         chromosome = chr_no,
                         from=chr_start, to=chr_end,
                         name = "treatment", 
                         col = "orange",
                         type='l')
ControlTrack <- DataTrack(range = control_list, 
                         chromosome = chr_no,
                         from=chr_start, to=chr_end,
                         name = "control", 
                         col = "blue",
                         type='l')
CTXTrack <- DataTrack(range = CTX_108, 
                         chromosome = chr_no,
                         from=chr_start, to=chr_end,
                         name = "Treatment-CTX108", 
                         col = "red",
                         type='h')
CTX1Track <- DataTrack(range = CTX_120, 
                      chromosome = chr_no,
                      from=chr_start, to=chr_end,
                      name = "control-CTX108", 
                      col = "black",
                      type='l')
plotTracks(c(CTXTrack,CTX1Track, 
             grTrack,
             g_axis),
           chromosome = chr_no, 
           from = chr_start, to = chr_end,
           extend.left=8000)
plotTracks(c(ControlTrack, 
             TreatmentTrack,
             grTrack,
             g_axis),
           chromosome = chr_no, 
           from = chr_start, to = chr_end,
           extend.left=8000)

# global_init.R
# Initializes global settings, logging, database connections, and color schemes.
# This script is sourced once at app launch to set up the environment for Shiny and plotting.

### Part0: Set up ###
requiredCRAN <- c('tidyverse','log4r')
requiredBiocPackages <- c('rtracklayer','GenomicFeatures','Gviz','ensembldb','Rsamtools','GenomicAlignments','glue')
purrr::walk(requiredCRAN, function(x) library(x,character.only = TRUE))
purrr::walk(requiredBiocPackages, function(x) library(x,character.only = TRUE))
my_layout <- function(level, ...) {
  paste(format(Sys.time()), "[", level, "]", ..., sep = " ", collapse = "", "\n")
}
logger <- logger(appenders=console_appender(my_layout))
level(logger) <- "INFO"

### Part1: Set Directory and corlor###
# data_dir <- "/mnt/gtklab01/xiaoqing/Run2_Nov_18/star/result"
# info(logger,glue("You have set your data directory as {data_dir}"))
treatment_group <- c(104,108,128,154)
control_group <- c(120,125,147,148)
info(logger,glue("You have defined CTX_{paste(treatment_group, collapse = ', CTX_')} as treatment and CTX_{paste(control_group, collapse = ', CTX_')} as control"))
condition_tag <- c("u","m1","m2","d")
control_shades <- generate_shades("#282A62")
treatment_shades <- generate_shades("#912C2C")

### Part2: Getting data ###
# Parsing GTF
# DB <- ensembldb::ensDbFromGtf("/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf",
#                               outfile="/mnt/gtklab01/xiaoqing/scaffold/GRCm39_Ens106.sqlite",
#                               organism = "Mus_musculus",
#                               genomeVersion ="GRCm39",
#                               version=106)
# save(DB,file="~/Capstone/test_data/db.rda")
load("~/Capstone/test_data/db.rda")
GRCm39 <- EnsDb(DB)
geneRanges_GRCm39 <- genes(GRCm39)
# Cytobands
if (!file.exists("/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz")) {
  info(logger,"The 'cytoBandIdeo.txt.gz' does not exist. Performing some action...")
  download.file(url="https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cytoBandIdeo.txt.gz",
              destfile="/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz")
} else {
  info(logger,"The 'cytoBandIdeo.txt.gz' has been downloaded.")
}
cytobands <- read_tsv("/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz",
         col_names = c("chrom","chromStart","chromEnd","name","gieStain"),show_col_types = FALSE)
# # get bigwig file path
# bw_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
#        group=c(control_group,treatment_group)),
#        tag=condition_tag) |>
#   mutate(bw_path=glue::glue("{data_dir}/unmapped_CTX_{group}_{tag}.bw"))
#get bam file path
# bam_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
#        group=c(control_group,treatment_group)),
#        tag=condition_tag) |>
#   mutate(bam_path=glue::glue("{data_dir}/unmappedAligned.sortedByCoord.out_CTX_{group}_{tag}.bam"))
#   "/mnt/gtklab01/xiaoqing/2025-01-14/analysis/bam/${sample%}_${tag}.bam"
#get bed file path
# intron_data <- import.bed("/mnt/gtklab01/xiaoqing/decoy/decoy3/intronic.bed")
# mashmap_data <- import.bed("/mnt/gtklab01/xiaoqing/decoy/decoy3/genome_found_sorted.bed")
# exon_data <- import.bed("/mnt/gtklab01/xiaoqing/decoy/decoy3/exon_out.bed")
options(ucscChromosomeNames=FALSE)
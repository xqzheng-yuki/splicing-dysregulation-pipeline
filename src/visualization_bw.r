# File: visualize_mapping.R
# Purpose: Visualize mapping results using GViz
# Dependencies: tidyverse, rtracklayer, GenomicFeatures, Gviz, ensembldb, Rsamtools, GenomicAlignments, BRGenomics, randomcoloR
# Inputs: bigwig files, GTF files
# Outputs: Visualization plots (PDF)
# Note: This code file is extracted from visualization_bw.qmd, please refer back for detailed information and explaination.

### Set up ###
requiredCRAN <- c('tidyverse','log4r')
requiredBiocPackages <- c('rtracklayer','GenomicFeatures','Gviz','ensembldb','Rsamtools','GenomicAlignments','BRGenomics','randomcoloR')
purrr::walk(requiredCRAN, function(x) library(x,character.only = TRUE))
purrr::walk(requiredBiocPackages, function(x) library(x,character.only = TRUE))
logger <- logger(appenders=console_appender(logfmt_log_layout()))
level(logger) <- "INFO"

### Parsing GTF ###
DB <- ensembldb::ensDbFromGtf("/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf",
                              outfile="/mnt/gtklab01/xiaoqing/scaffold/GRCm39_Ens106.sqlite",
                              organism = "Mus_musculus",
                              genomeVersion ="GRCm39",
                              version=106)
GRCm39 <- EnsDb(DB)
# seqlevelsStyle(ens106) <- "UCSC"
# GRCm39 <- ens106
geneRanges_GRCm39 <- genes(GRCm39)

### Cytobands ###
if (!file.exists("/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz")) {
  info(logger,"The 'cytoBandIdeo.txt.gz' does not exist. Performing some action...")
  download.file(url="https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cytoBandIdeo.txt.gz",
              destfile="/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz")
} else {
  info(logger,"The 'cytoBandIdeo.txt.gz' has been downloaded.")
}
cytobands <- read_tsv("/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz",
         col_names = c("chrom","chromStart","chromEnd","name","gieStain"),show_col_types = FALSE)

### Take in bigwig file ###
# self-define process function #
source("~/Capstone/src/dependent_bw.r")

# parameter #
data_dir <- "/mnt/gtklab01/xiaoqing/star/results/group/Nov_18"

treatment_group <- c(104,108,128,154)
control_group <- c(120,125,147,148)
condition_tag <- c("u","m1","m2","d")

control_shades <- generate_shades("#282A62")
treatment_shades <- generate_shades("#912C2C")

# get bigwig file path #
bw_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
       group=c(control_group,treatment_group)),
       tag=condition_tag) |>
  mutate(bw_path=glue::glue("{data_dir}/unmapped_CTX_{group}_{tag}.bw"))

#get bam file path #
bam_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
       group=c(control_group,treatment_group)),
       tag=condition_tag) |>
  mutate(bam_path=glue::glue("{data_dir}/unmappedAligned.sortedByCoord.out_CTX_{group}_{tag}.bam"))

#get bed file path #
intron_data <- import.bed("/mnt/gtklab01/xiaoqing/salmon/index/decoy3/intronic.bed")
mashmap_data <- import.bed("/mnt/gtklab01/xiaoqing/salmon/index/decoy3/genome_found_sorted.bed")
options(ucscChromosomeNames=FALSE)
#########################################################
################### plotting function ###################
#########################################################
plot_gene <- function(goi) {
  
  if (!startsWith(goi,"ENSMUSG")) {
  debug(logger,"Not inputing gene id, trying to change into gene id")
  goi <- get_gene_id(goi)
  }

  gr <- geneRanges_GRCm39[goi]
  
  ## let's make it 10% wider
  embiggen_factor <- round(width(gr)*0.1)
  start(gr) <- start(gr) - embiggen_factor
  end(gr) <- end(gr) + embiggen_factor
 
  ## first track:: cytobands
  idt <- IdeogramTrack(chromosome=as.character(chrom(gr)),
                     genome='mm39',
                     bands=cytobands,
                     from=start(gr),to=end(gr))
  show_track <- add_track(idt)

  ## track:: current trxns
  grTrack <- GeneRegionTrack(ensembldb::filter(GRCm39, ~ gene_id == goi),
                             name=gr$symbol,
                             start=start(gr),
                             end=end(gr),
                             cex.group=0.4,
                             showId=TRUE)
  show_track <- add_track(grTrack,show_track)
  
  ## track:: for bed
  intron_gr <- intron_data[queryHits(findOverlaps(intron_data, gr))]
  intron_track <- AnnotationTrack(intron_gr,
                  chromosome=as.character(chrom(gr)),
                  genome='mm39',
                  name = "IN",
                  stacking = "dense",
                  start=start(gr),end=end(gr))
  show_track <- add_track(intron_track,show_track)
  mashmap_gr <- mashmap_data[queryHits(findOverlaps(mashmap_data, gr))]
  mashmap_track <- AnnotationTrack(mashmap_gr,
                  chromosome=as.character(chrom(gr)),
                  genome='mm39',
                  name = "MM",
                  stacking = "dense",
                  start=start(gr),end=end(gr))
  show_track <- add_track(mashmap_track,show_track)
  debug(logger,sprintf("length of show track is %d",length(show_track)))
  ## track:: for bw
  debug(logger,sprintf("length of gr is %d",length(gr)))
  control_bw_data <- merge_treatment_bw(selection_range = gr,
                                     treatment_n = "control")
  ko_bw_data <- merge_treatment_bw(selection_range = gr,
                                     treatment_n = "treatment")
  ylim <- get_ylim(first_set=control_bw_data,second_set=ko_bw_data)
  control_datatrack <- map(names(control_bw_data),
    function(n) {
      color <- control_shades[match(n,names(control_bw_data))]
      DataTrack(control_bw_data[[n]],
                name = glue::glue("control_{n}"),
                type="polygon",
                ylim=ylim,
                fill.mountain=c("white",color),
                col="black")
    })
  show_track <- add_track(control_datatrack,show_track)

  ko_datatrack <- map(names(ko_bw_data),
    function(n) {
      color <- treatment_shades[match(n,names(ko_bw_data))]
      DataTrack(ko_bw_data[[n]],
                name = glue::glue("treatment_{n}"),
                type="polygon",
                ylim=ylim,
                fill.mountain=c("white",color),
                col="black",stackHeight=1)
    })
  show_track <- add_track(ko_datatrack,show_track)
  
  ## track:: genomic coordinates
  show_track <- add_track(GenomeAxisTrack(),show_track)
  pdf(file = paste0("~/Capstone/results/", gr$symbol, ".pdf"))
  plotTracks(show_track, cex.sampleNames = 0.5, main = gr$symbol, fontface.main = 1.5)
  dev.off()
}


plot_supplymental <- function(goi) {

  if (!startsWith(goi,"ENSMUSG")) {
  debug(logger,"Not inputing gene id, trying to change into gene id")
  goi <- get_gene_id(goi)
  }

  gr <- geneRanges_GRCm39[goi]
  
  ## let's make it 10% wider
  embiggen_factor <- round(width(gr)*0.1)
  start(gr) <- start(gr) - embiggen_factor
  end(gr) <- end(gr) + embiggen_factor

  ## track:: current trxns
  grTrack <- GeneRegionTrack(ensembldb::filter(GRCm39, ~ gene_id == goi),
                             name=gr$symbol,
                             start=start(gr),
                             end=end(gr),
                             cex.group=0.4,
                             showId=TRUE)
  show_track <- add_track(grTrack)
  
  ## track:: for bw
  control_bw_data <- merge_treatment_bw(selection_range = gr,
                                     treatment_n = "control")
  ko_bw_data <- merge_treatment_bw(selection_range = gr,
                                     treatment_n = "treatment")
  ylim_list <- lapply(names(control_bw_data), function(condition_tag) {
    get_ylim(control_bw_data, ko_bw_data, condition_tag)
  })
  control_datatrack <- map(names(control_bw_data),
    function(n) {
      color <- control_shades[match(n,names(control_bw_data))]
      ylim <- ylim_list[match(n,names(control_bw_data))]
      DataTrack(control_bw_data[[n]],
                name = glue::glue("control_{n}"),
                type="polygon",
                fill.mountain=c("white",color),
                col="black")
    })

  ko_datatrack <- map(names(ko_bw_data),
    function(n) {
      color <- treatment_shades[match(n,names(ko_bw_data))]
      ylim <- ylim_list[match(n,names(ko_bw_data))]
      DataTrack(ko_bw_data[[n]],
                name = glue::glue("treatment_{n}"),
                type="polygon",
                fill.mountain=c("white",color),
                col="black")
    })

  
  stack_datatrack <- map(names(control_bw_data),
    function(n) {
      index <- match(n,names(control_bw_data))
      ylim <- ylim_list[index]
      title <- paste0("stack_",n)
      debug(logger,title)
      OverlayTrack(trackList = list(control_datatrack[[index]],
                                    ko_datatrack[[index]]),
                  name = title,
                  ylim = ylim
      )
    })
  show_track <- add_track(stack_datatrack,show_track)
  
  alignment_track <- pmap(bam_fileinfo, function(treatment,group,tag,bam_path) {
    AlignmentsTrack(bam_path, start = start(gr), end= end(gr),chromosome=as.character(chrom(gr)),
    type=c("coverage", "sashimi"),name=paste0("CTX_",group,"_",tag))
  })
  m2row <- which(bam_fileinfo$tag=="m2")

  ## track:: genomic coordinates
  show_track <- add_track(GenomeAxisTrack(),show_track)

  options(ucscChromosomeNames=FALSE)
  pdf(file = paste0("~/Capstone/results/", gr$symbol, "_supplymental.pdf"))
  plotTracks(show_track, cex.sampleNames = 0.5, main = gr$symbol, fontface.main = 1.5)
  plotTracks(alignment_track[m2row], main = "m2",from = start(gr), to = end(gr))
  dev.off()
}
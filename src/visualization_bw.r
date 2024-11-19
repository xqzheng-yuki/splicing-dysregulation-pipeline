# File: visualize_mapping.R
# Purpose: Visualize mapping results using GViz
# Dependencies: tidyverse, rtracklayer, GenomicFeatures, Gviz, ensembldb, Rsamtools, GenomicAlignments, BRGenomics, randomcoloR
# Inputs: bigwig files, GTF files
# Outputs: Visualization plots (PDF)
# Note: This code file is extracted from visualization_bw.qmd, please refer back for detailed information and explaination.

### Set up ###
requiredCRAN <- c('tidyverse')
requiredBiocPackages <- c('rtracklayer','GenomicFeatures','Gviz','ensembldb','Rsamtools','GenomicAlignments','BRGenomics','randomcoloR')
purrr::walk(requiredCRAN, function(x) library(x,character.only = TRUE))
purrr::walk(requiredBiocPackages, function(x) library(x,character.only = TRUE))

### Parsing GTF ###
DB <- ensembldb::ensDbFromGtf("/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf",
                              outfile="/mnt/gtklab01/xiaoqing/scaffold/GRCm39_Ens106.sqlite",
                              organism = "Mus_musculus",
                              genomeVersion ="GRCm39",
                              version=106)
ens106 <- EnsDb(DB)
seqlevelsStyle(ens106) <- "UCSC"
GRCm39 <- ens106
geneRanges_GRCm39 <- genes(GRCm39)

### Cytobands ###
if (!file.exists("/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz")) {
  message("The 'cytoBandIdeo.txt.gz' does not exist. Performing some action...")
  download.file(url="https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cytoBandIdeo.txt.gz",
              destfile="/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz")
} else {
  message("The 'cytoBandIdeo.txt.gz' has been downloaded.")
}
cytobands <- read_tsv("/mnt/gtklab01/xiaoqing/scaffold/cytoBandIdeo.txt.gz",
         col_names = c("chrom","chromStart","chromEnd","name","gieStain"))

### Take in bigwig file ###
# self-define process function #
import_bigwig <- function(file,bw_selection) {
  bw_data <- BRGenomics::makeGRangesBRG(import.bw(file,selection=rtracklayer::BigWigSelection(bw_selection)))
  keepStandardChromosomes(bw_data,pruning.mode = "coarse")
}
merge_treatment_bw <- function(selection_range,treatment_n,bwinfo=bw_fileinfo) {
  bwinfo <- dplyr::filter(bwinfo,treatment==treatment_n)
  bw_data <- split(map(bwinfo$bw_path,import_bigwig,bw_selection=selection_range),
        bwinfo$tag)
  bw_range <- map(bw_data, ~ {
                  mergeGRangesData(.x,exact_overlaps = TRUE)
    })
  return(bw_range)
}
generate_shades <- function(base_color,condition_tag) {
  # Convert hex to RGB
  base_rgb <- col2rgb(base_color) / 255
  shades <- colorRampPalette(c("white", rgb(base_rgb[1], base_rgb[2], base_rgb[3])))(5)
  return(shades[2:length(shades)])
}
add_track <- function(new_track,tracks=list()) {
  if (is.list(new_track)) {
    tracks <- append(tracks,new_track)
  } else {
  tracks <- append(tracks,list(new_track))
  }
  return(tracks)
}
get_gene_id <- function(geneName) {
  gene_id_match <- geneRanges_GRCm39[which(elementMetadata(geneRanges_GRCm39)$gene_name==geneName)]$gene_id
  if (length(gene_id_match) == 0) {
    message(paste0("No matches found for gene name ",geneName,"."))
    message("Testing for approximate results...")
    approx_name <- str_to_title(geneName)
    gene_id_match <- geneRanges_GRCm39[which(elementMetadata(geneRanges_GRCm39)$gene_name==approx_name)]$gene_id
    if (length(gene_id_match) != 0) {
      message(paste0("Found approximate match with gene name ",approx_name))
    } else {
      message("Not found approximately matched gene name. Please double check the gene name or provide a gene id starting with ENSMUSG.")
      return(invisible())
    }
  }
  message(paste0("The gene id is ",gene_id_match,"."))
  return(gene_id_match)
}
get_gene_id("Ppp6c")
get_gene_id("UNC13A")
get_gene_id("UNC13")
# parameter #
data_dir <- "/mnt/gtklab01/xiaoqing/star/results/group/Nov_18"

treatment_group <- c(104,108,128,154)
control_group <- c(120,125,147,148)
condition_tag <- c("u","m1","m2","d")

control_shades <- generate_shades("#282A62")
treatment_shades <- generate_shades("#912C2C")
# get bigwig file information #
bw_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
       group=c(control_group,treatment_group)),
       tag=condition_tag) |>
  mutate(bw_path=glue::glue("{data_dir}/unmapped_CTX_{group}_{tag}.bw"))
# plotting function #
plot_gene <- function(goi) {
  
  if (!startsWith(goi,"ENSMUSG")) {
  message("Not inputing gene id, trying to change into gene id")
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
  
  ## track:: for bw
  control_bw_data <-merge_treatment_bw(selection_range = gr,
                                     treatment_n = "control")
  ko_bw_data <-merge_treatment_bw(selection_range = gr,
                                     treatment_n = "treatment")

  control_datatrack <- map(names(control_bw_data),
    function(n) {
      color <- control_shades[match(n,names(control_bw_data))]
      DataTrack(control_bw_data[[n]],
                name = glue::glue("control_{n}"),
                type="polygon",
                # ylim=c(0,600),
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
                # ylim=c(0,600),
                fill.mountain=c("white",color),
                col="black")
    })
  show_track <- add_track(ko_datatrack,show_track)
  
  ## track:: genomic coordinates
  show_track <- add_track(GenomeAxisTrack(),show_track)
  pdf(file = paste0("~/Capstone/results/", gr$symbol, ".pdf"))
  plotTracks(show_track, cex.sampleNames = 0.6, main = gr$symbol)
  dev.off()
}

plot_gene("ENSMUSG00000026753.6")
plot_gene("Ube2d1")
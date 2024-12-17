# File: visualize_mapping.R
# Purpose: Visualize mapping shiny component 
# Dependencies: tidyverse, rtracklayer, GenomicFeatures, Gviz, ensembldb, Rsamtools, GenomicAlignments, BRGenomics, randomcoloR
# Inputs: bigwig files, GTF files
# Note: This code file is extracted from visualization_bw.qmd, and slightly modified for shiny app.
# please refer back for detailed information and explaination.



### Part0: Set up ###
requiredCRAN <- c('tidyverse','log4r')
requiredBiocPackages <- c('rtracklayer','GenomicFeatures','Gviz','ensembldb','Rsamtools','GenomicAlignments','BRGenomics','randomcoloR','glue')
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
# #get bam file path
# bam_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
#        group=c(control_group,treatment_group)),
#        tag=condition_tag) |>
#   mutate(bam_path=glue::glue("{data_dir}/unmappedAligned.sortedByCoord.out_CTX_{group}_{tag}.bam"))
#get bed file path
# intron_data <- import.bed("/mnt/gtklab01/xiaoqing/decoy/decoy3/intronic.bed")
# mashmap_data <- import.bed("/mnt/gtklab01/xiaoqing/decoy/decoy3/genome_found_sorted.bed")
# exon_data <- import.bed("/mnt/gtklab01/xiaoqing/decoy/decoy3/exon_out.bed")
options(ucscChromosomeNames=FALSE)

### Part3: Getting Track Information ###
tracklist <- function(goi,run_number){
  debug(logger,glue("You're now working on extracting tracks for visualization."))
  run_database <- get_data_dir(run_number)
  run_decoy <- get_decoy_dir(run_number)
  bw_fileinfo <- get_bw_path(run_database)
  bam_fileinfo <- get_bam_path(run_database)
  intron_data <- get_intron_bed(run_decoy)
  mashmap_data <- get_mashmap_bed(run_decoy)
  exon_data <- get_exon_bed(run_decoy)
  gr <- geneRanges_GRCm39[goi]

    ## track:: current trxns
    grTrack <- GeneRegionTrack(ensembldb::filter(GRCm39, ~ gene_id == goi),
                             name=gr$symbol,
                             start=start(gr),
                             end=end(gr),
                             cex.group=0.4,
                             showId=TRUE)
    show_track <- add_track(grTrack)
    debug(logger,glue("Added transcriptome track."))
    debug(logger,glue("Length of show track now is {length(show_track)}"))
  
    ## track:: for bed
    exon_gr <- exon_data[queryHits(findOverlaps(exon_data, gr))]
    exon_track <- AnnotationTrack(exon_gr,
                                  chromosome=as.character(chrom(gr)),
                                  genome='mm39',
                                  name = "MM",
                                  stacking = "dense",
                                  start=start(gr),end=end(gr),
                                  fill='blue')
    show_track <- add_track(exon_track,show_track)
    debug(logger,glue("Added exon region track. There are {length(exon_track)} exon(s)."))
    debug(logger,glue("Length of show track now is {length(show_track)}"))
    intron_gr <- intron_data[queryHits(findOverlaps(intron_data, gr))]
    intron_track <- AnnotationTrack(intron_gr,
                  chromosome=as.character(chrom(gr)),
                  genome='mm39',
                  name = "IN",
                  stacking = "dense",
                  start=start(gr),end=end(gr),
                  fill='grey')
    show_track <- add_track(intron_track,show_track)
    debug(logger,glue("Added pure intronic region track. There are {length(intron_track)} intronic(s)."))
    debug(logger,glue("Length of show track now is {length(show_track)}"))
    seqlevels(mashmap_data, pruning.mode = "coarse") <- seqlevels(gr)
    mashmap_gr <- mashmap_data[queryHits(findOverlaps(mashmap_data, gr))]
    mashmap_track <- AnnotationTrack(mashmap_gr,
                  chromosome=as.character(chrom(gr)),
                  genome='mm39',
                  name = "MM",
                  stacking = "dense",
                  start=start(gr),end=end(gr),
                  fill='#FFA07A')
    show_track <- add_track(mashmap_track,show_track)
    debug(logger,glue("Added mashmap alignment region"))
    debug(logger,glue("Length of show track now is {length(show_track)}"))

    ## track:: for bw
    debug(logger,glue("There are {length(rownames(bw_fileinfo))} lines of 'bw_fileinfo', right before 'control_bw_data'."))
    control_bw_data <- merge_treatment_bw(selection_range = gr,
                                     treatment_n = "control",
                                     bwinfo = bw_fileinfo)
    debug(logger,glue("Now the bw_fileinfo is {length(rownames(bw_fileinfo))}, right before 'ko_bw_data'."))
    ko_bw_data <- merge_treatment_bw(selection_range = gr,
                                     treatment_n = "treatment",
                                     bwinfo = bw_fileinfo)
    ylim <- get_ylim(first_set=control_bw_data,second_set=ko_bw_data)
    ylim_list <- lapply(names(control_bw_data), function(condition_tag) {
        get_ylim(control_bw_data, ko_bw_data, condition_tag)
    })
    debug(logger,glue("There are {length(rownames(bw_fileinfo))} lines of 'bw_fileinfo', right before 'control_datatrack'."))
    control_datatrack <- map(names(control_bw_data),
        function(n) {
            color <- control_shades[match(n,names(control_bw_data))]
            DataTrack(control_bw_data[[n]],
                name = glue::glue("ctrl_{n}"),
                type="polygon",
                ylim = ylim,
                fill.mountain=c("white",color),
                col="black")
        })
    show_track <- add_track(control_datatrack,show_track)
    debug(logger,glue("Added {length(control_datatrack)} track(s) bw data for control group"))
    debug(logger,glue("Length of show track now is {length(show_track)}"))
    debug(logger,glue("There are {length(rownames(bw_fileinfo))} lines of 'bw_fileinfo', right before 'ko_datatrack'."))
    ko_datatrack <- map(names(ko_bw_data),
        function(n) {
            color <- treatment_shades[match(n,names(ko_bw_data))]
            DataTrack(ko_bw_data[[n]],
                name = glue::glue("ko_{n}"),
                type="polygon",
                ylim = ylim,
                fill.mountain=c("white",color),
                col="black")
        })
    show_track <- add_track(ko_datatrack,show_track)
    debug(logger,glue("Added {length(ko_datatrack)} track(s) bw data for treatment group"))
    debug(logger,glue("Length of show track now is {length(show_track)}"))
    
    ## track:: for overlap track
    stack_datatrack <- map(names(control_bw_data),
        function(n) {
            index <- match(n,names(control_bw_data))
            ylim <- ylim_list[index]
            title <- paste0("stack_",n)
            OverlayTrack(trackList = list(control_datatrack[[index]],ko_datatrack[[index]]),
                name = title,
                ylim = ylim
            )
    })
    show_track <- add_track(stack_datatrack,show_track)
    debug(logger,glue("Added {length(stack_datatrack)} track(s) overlay bw data grouped by type of unmapped"))
    debug(logger,glue("Length of show track now is {length(show_track)}"))

    ## track:: sashimi
    alignment_track <- pmap(bam_fileinfo, function(treatment,group,tag,bam_path) {
        AlignmentsTrack(bam_path, start = start(gr), end= end(gr),chromosome=as.character(chrom(gr)),
        type=c("coverage", "sashimi"),name=paste0(group,"_",tag))
    })
    drow <- which(bam_fileinfo$tag=="d")
    m1row <- which(bam_fileinfo$tag=="m1")
    m2row <- which(bam_fileinfo$tag=="m2")
    d_ali_track <- alignment_track[drow]
    m1_ali_track <- alignment_track[m1row]
    m2_ali_track <- alignment_track[m2row]
    show_track <- add_track(d_ali_track,show_track)
    debug(logger,glue("Added {length(d_ali_track)} track(s) sashimi plot for map to decoy"))
    debug(logger,glue("Length of show track now is {length(show_track)}"))
    show_track <- add_track(m1_ali_track,show_track)
    debug(logger,glue("Added {length(m1_ali_track)} track(s) sashimi plot for map to m1"))
    debug(logger,glue("Length of show track now is {length(show_track)}"))
    show_track <- add_track(m2_ali_track,show_track)
    debug(logger,glue("Added {length(m2_ali_track)} track(s) sashimi plot for map to m2"))
    debug(logger,glue("Length of show track now is {length(show_track)}"))

    ## track:: genomic coordinates
    show_track <- add_track(GenomeAxisTrack(),show_track)
    debug(logger,glue("Added 1 track(s) genomic coordinate axis"))
    info(logger,glue("Length of show track now is {length(show_track)}"))
    return(show_track)
}

plotplot <- function(tracklist,goi) {
    gr <- enlarge_gr(goi)
    plotTracks(tracklist, cex.sampleNames = 0.8, from = start(gr), to = end(gr), main = gr$symbol, fontface.main = 1.5, fontsize = 15)
}

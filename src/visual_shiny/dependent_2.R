source("~/Capstone/src/visual_shiny/dependent_bw.r")
source("~/Capstone/src/visual_shiny/plot_shiny.r")
source("~/Capstone/src/visual_shiny/major_plot_function.r")
filter_bam_path <- function(path) {
  data_dir <- path
  bam_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                     group=c(control_group,treatment_group)),tag="d") |>
    mutate(bam_path=glue::glue("{data_dir}/CTX_{group}_d.bam"))
  return(bam_fileinfo)
}
filter_bw_path <- function(path) {
  data_dir <- path
  bw_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                     group=c(control_group,treatment_group)),tag="d") |>
    mutate(bw_path=glue::glue("{data_dir}/bw/CTX_{group}_d.bw"))
  return(bw_fileinfo)
}

get_gene_name <- function(geneID){
  matches <- which(elementMetadata(geneRanges_GRCm39)$gene_id == geneID)
  if (length(matches) == 0) {
    #  warn(logger(),"There is no match found.")
    gene_name_match = NA
  } else {
    gene_name_match <- geneRanges_GRCm39[matches]$gene_name
  }
  # info(logger,paste0("The gene ",geneID," match is ",gene_name_match))
  return(gene_name_match)
}

trackset <- function(goi){
  debug(logger,glue("You're now working on extracting tracks for visualization(II)."))
  run_database <- "/mnt/gtklab01/xiaoqing/2025-01-14/filter"
  run_decoy <- "/mnt/gtklab01/xiaoqing"
  # bw_fileinfo <- filter_bw_path(run_database)
  bam_fileinfo <- filter_bam_path(run_database)
  # intron_data <- get_intron_bed(run_decoy)
  exon_data <- get_exon_bed(run_decoy)
  gr <- geneRanges_GRCm39[goi]
  
  ## track:: current trxns
  grTrack <-  GeneRegionTrack(ensembldb::filter(GRCm39, ~ gene_id == goi),
                             name=gr$symbol,
                             start=start(gr),
                             end=end(gr),
                             cex.group=0.4,
                             showId=TRUE,
                             stacking="dense")
  show_track <- add_track(grTrack)
  debug(logger,glue("Added transcriptome track."))
  debug(logger,glue("Length of show track now is {length(show_track)}"))
  
  ## track:: for bed
  exon_gr <- exon_data[queryHits(findOverlaps(exon_data, gr))]
  exon_track <- AnnotationTrack(exon_gr,
                                chromosome=as.character(chrom(gr)),
                                genome='mm39',
                                name = "exon",
                                stacking = "dense",
                                start=start(gr),end=end(gr),
                                fill='blue')
  show_track <- add_track(exon_track,show_track)
  debug(logger,glue("Added exon region track. There are {length(exon_track)} exon(s)."))
  debug(logger,glue("Length of show track now is {length(show_track)}"))
  
  ## track:: for bw
  # control_bw_data <- merge_treatment_bw(selection_range = gr,
  #                                       treatment_n = "control",
  #                                       bwinfo = bw_fileinfo)
  # ko_bw_data <- merge_treatment_bw(selection_range = gr,
  #                                  treatment_n = "treatment",
  #                                  bwinfo = bw_fileinfo)
  # ylim <- get_ylim(first_set=control_bw_data,second_set=ko_bw_data)
  # control_datatrack <- map(names(control_bw_data),
  #                          function(n) {
  #                            color <- control_shades[match(n,names(control_bw_data))]
  #                            DataTrack(control_bw_data[[n]],
  #                                      name = glue::glue("ctrl_{n}"),
  #                                      type="polygon",
  #                                      ylim = ylim,
  #                                      fill.mountain=c("white",color),
  #                                      col="black")
  #                          })
  # show_track <- add_track(control_datatrack,show_track)
  # ko_datatrack <- map(names(ko_bw_data),
  #                     function(n) {
  #                       color <- treatment_shades[match(n,names(ko_bw_data))]
  #                       DataTrack(ko_bw_data[[n]],
  #                                 name = glue::glue("ko_{n}"),
  #                                 type="polygon",
  #                                 ylim = ylim,
  #                                 fill.mountain=c("white",color),
  #                                 col="black")
  #                     })
  # show_track <- add_track(ko_datatrack,show_track)

  ## track:: sashimi
  alignment_track <- pmap(bam_fileinfo, function(treatment,group,tag,bam_path) {
    AlignmentsTrack(bam_path, start = start(gr), end= end(gr),chromosome=as.character(chrom(gr)),
                    type=c("coverage", "sashimi"),name=paste0("CTX_",group))
  })
  show_track <- add_track(alignment_track,show_track)
  debug(logger,glue("Added {length(alignment_track)} track(s) sashimi plot for decoy pure"))
  debug(logger,glue("Length of show track now is {length(show_track)}"))
  
  ## track:: genomic coordinates
  show_track <- add_track(GenomeAxisTrack(),show_track)
  debug(logger,glue("Added 1 track(s) genomic coordinate axis"))
  info(logger,glue("Length of show track now is {length(show_track)}"))
  return(show_track)
}

dev.new(width = 13, height = 30,noRStudioGD = T,file="/mnt/gtklab01/xiaoqing/analysis/graph/Only_in_spliceat.pdf")
plotplot(trackset(get_gene_id("Evc2")),get_gene_id("Evc2"))
plotplot(trackset(get_gene_id("Rps13")),get_gene_id("Rps13"))
plotplot(trackset(get_gene_id("Ccdc187")),get_gene_id("Ccdc187"))
dev.off()

dev.new(width = 13, height = 30,noRStudioGD = T,file="/mnt/gtklab01/xiaoqing/analysis/graph/Only_in_decap.pdf")
plotplot(trackset(get_gene_id("Ptp4a2")),get_gene_id("Ptp4a2"))
plotplot(trackset(get_gene_id("Mfsd13a")),get_gene_id("Mfsd13a"))
plotplot(trackset(get_gene_id("Heatr5b")),get_gene_id("Heatr5b"))
dev.off()
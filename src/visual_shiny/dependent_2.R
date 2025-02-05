filter_bam_path <- function(path) {
  data_dir <- path
  bam_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                     group=c(control_group,treatment_group))) |>
    mutate(bam_path=glue::glue("{data_dir}/CTX_{group}_d.bam"))
  return(bam_fileinfo)
}

trackset <- function(goi){
  debug(logger,glue("You're now working on extracting tracks for visualization(II)."))
  run_database <- "/mnt/gtklab01/xiaoqing/2025-01-14/filter"
  run_decoy <- "/mnt/gtklab01/xiaoqing/decoys/2025-01-14"
  bam_fileinfo <- filter_bam_path(run_database)
  # intron_data <- get_intron_bed(run_decoy)
  exon_data <- get_exon_bed(run_decoy)
  gr <- geneRanges_GRCm39[goi]
  
  ## track:: current trxns
  grTrack <- GeneRegionTrack(ensembldb::filter(GRCm39, ~ gene_id == goi),
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
                                name = "MM",
                                stacking = "dense",
                                start=start(gr),end=end(gr),
                                fill='blue')
  show_track <- add_track(exon_track,show_track)
  debug(logger,glue("Added exon region track. There are {length(exon_track)} exon(s)."))
  debug(logger,glue("Length of show track now is {length(show_track)}"))
  # intron_gr <- intron_data[queryHits(findOverlaps(intron_data, gr))]
  # intron_track <- AnnotationTrack(intron_gr,
  #                                 chromosome=as.character(chrom(gr)),
  #                                 genome='mm39',
  #                                 name = "IN",
  #                                 stacking = "dense",
  #                                 start=start(gr),end=end(gr),
  #                                 fill='grey')
  # show_track <- add_track(intron_track,show_track)
  # debug(logger,glue("Added pure intronic region track. There are {length(intron_track)} intronic(s)."))
  # debug(logger,glue("Length of show track now is {length(show_track)}"))
  # 
  ## track:: sashimi
  alignment_track <- pmap(bam_fileinfo, function(treatment,group,bam_path) {
    AlignmentsTrack(bam_path, start = start(gr), end= end(gr),chromosome=as.character(chrom(gr)),
                    type=c("sashimi"),name=paste0(group))
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

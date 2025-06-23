source("~/Capstone/src/visual_shiny/dependent_bw.r")
source("~/Capstone/src/visual_shiny/plot_shiny.r")
source("~/Capstone/src/visual_shiny/major_plot_function.r")
filter_bam_path <- function(path) {
  data_dir <- path
  bam_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                     group=c(control_group,treatment_group)),tag="d") |>
    # mutate(bam_path=glue::glue("{data_dir}/CTX_{group}_d.bam"))
    mutate(bam_path=glue::glue("{data_dir}/CTX_{group}/d.sorted.bam"))
  return(bam_fileinfo)
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
  run_database <- "/mnt/gtklab01/xiaoqing/star/result"
  # run_database <- "/mnt/gtklab01/xiaoqing/analysis/bam"
  run_decoy <- "/mnt/gtklab01/xiaoqing"
  bam_fileinfo <- filter_bam_path(run_database)
  # run_database <- "/mnt/gtklab01/xiaoqing/Run4_Dec_09/star/result"
  # run_decoy <- "/mnt/gtklab01/xiaoqing/Z-decoy/decoy5"
  # bam_fileinfo <- get_bam_path(run_database)
  intron_data <- get_intron_bed(run_decoy)
  exon_data <- get_exon_bed(run_decoy)
  gr <- geneRanges_GRCm39[goi]
  
  ## track:: current trxns
  grTrack <-  GeneRegionTrack(ensembldb::filter(GRCm39, ~ gene_id == goi),
                             name=gr$symbol,
                             start=start(gr),
                             end=end(gr),
                             cex.group=0.4,
                             showId=TRUE,
                             stacking="full")
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

  ## track:: sashimi
#   alignment_track <- pmap(bam_fileinfo, function(treatment, group, tag, bam_path) {
#     tryCatch({
#       tryCatch({
#   at <- AlignmentsTrack(
#     bam_path,
#     start = start(gr), end = end(gr),
#     chromosome = as.character(chrom(gr)),
#     name = paste0("CTX_", group)
#   )
#   
#   # Check if it has any junctions (ranges) before enabling sashimi
#   if (length(ranges(at)) > 0) {
#     displayPars(at) <- list(type = c("coverage", "sashimi"))
#   } else {
#     message(glue("⚠️ Track {group} has no junctions — using coverage only"))
#     displayPars(at) <- list(type = "coverage")
#   }
#   
#   at
# }, error = function(e) {
#   warning(glue("⚠️ Failed to create AlignmentsTrack for {group}: {e$message}"))
#   NULL
# })
#     }, error = function(e) {
#       warning(logger, glue("⚠️ Failed to create sashimi track for {group}: {e$message}"))
#       NULL
#     })
#   })
#   alignment_track <- alignment_track[!sapply(alignment_track, is.null)]
#   if (length(alignment_track) > 0) {
#     show_track <- add_track(alignment_track, show_track)
#     debug(logger, glue("Added {length(alignment_track)} sashimi track(s)"))
#   } else {
#     warning(glue("⚠️ No valid sashimi tracks created for {goi}"))
#   }
  alignment_track <- pmap(bam_fileinfo, function(treatment,group,tag,bam_path) {
    AlignmentsTrack(bam_path, start = start(gr), end= end(gr),chromosome=as.character(chrom(gr)),
                    type=c("coverage","sashimi"),name=paste0(group,"_",tag))
  })
  show_track <- add_track(alignment_track,show_track)
  debug(logger, glue("Length of show track now is {length(show_track)}"))
  
  ## track:: genomic coordinates
  show_track <- add_track(GenomeAxisTrack(), show_track)
  debug(logger, glue("Added 1 track(s) genomic coordinate axis"))
  info(logger, glue("Length of show track now is {length(show_track)}"))

  info(logger, glue("ALL TRACKS AVAILABLE NOW!"))
  return(show_track)
}

  debug(logger, glue("Length of show track now is {length(show_track)}"))
  
  ## track:: genomic coordinates
  show_track <- add_track(GenomeAxisTrack(), show_track)
  debug(logger, glue("Added 1 track(s) genomic coordinate axis"))
  info(logger, glue("Length of show track now is {length(show_track)}"))

  info(logger, glue("ALL TRACKS AVAILABLE NOW!"))
  return(show_track)
}

save_result <- function(one){
  path <- paste0("/mnt/gtklab01/xiaoqing/result_for_thesis/SALMONcloseup_",one,".png")
  png(path, width = 500, height = 780, units = "px", pointsize = 12)
  # pdf(path, width = 13, height = 20)
  plotplot(trackset(get_gene_id(one)),get_gene_id(one))
  dev.off()
  gc()
}
  
only_in_splice <- c("Evc2","Rps13","Ccdc187")
only_in_decap <- c("Ptp4a2","Mfsd13a","Heatr5b")
common_one <-c("Ppp6c","Trim8","1110017D15Rik")
for (name in c(common_one,only_in_decap,only_in_splice)) {
  print(name)
  save_result(name)
}

# pdf("/mnt/gtklab01/xiaoqing/result_for_thesis/star_Only_in_spliceat.pdf", width = 13, height = 30)
plotplot(trackset(get_gene_id("Evc2")),get_gene_id("Evc2"))
# gc()
plotplot(trackset(get_gene_id("Rps13")),get_gene_id("Rps13"))
# gc()
plotplot(trackset(get_gene_id("Ccdc187")),get_gene_id("Ccdc187"))
# gc()
# dev.off()
# 
# pdf("/mnt/gtklab01/xiaoqing/result_for_thesis/star_Only_in_decap.pdf", width = 13, height = 30)
plotplot(trackset(get_gene_id("Ptp4a2")),get_gene_id("Ptp4a2"))
# gc()
plotplot(trackset(get_gene_id("Mfsd13a")),get_gene_id("Mfsd13a"))
# gc()
plotplot(trackset(get_gene_id("Heatr5b")),get_gene_id("Heatr5b"))
# gc()
# dev.off()
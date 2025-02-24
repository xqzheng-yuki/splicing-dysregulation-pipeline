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
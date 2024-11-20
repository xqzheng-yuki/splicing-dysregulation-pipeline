selected_bam_files <- function(treatment,tag,fileinfo=bam_fileinfo) {
    .tag = tag
    .treatment = treatment
    dplyr::filter(fileinfo,tag==.tag,treatment==.treatment) |> pull(bam_path)
}
treated_m1 <- selected_bam_files('treatment','m1')

import_bam <- function(file,bam_selection) {
  bam_selection <- gr
  bam_param <- ScanBamParam(which = bam_selection)
  bam_data <- GenomicAlignments::readGAlignmentPairs(file, param = bam_param)
  genome(bam_data) <- genome(GRCm39)
  bam_data
}

combined_m1 <- do.call("c" ,purrr::map(treated_m1,import_bam,gr))
combined_m1

export(combined_m1,con="combined_m1.bam")

test_combined_track <- AlignmentsTrack("combined_m1.bam",start = start(gr), end= end(gr),chromosome="chr2",
    type=c("coverage", "sashimi"),name=paste0("CTX_"))
plotTracks(test_combined_track, from = start(gr), to = end(gr),chrom="chr2")

bam_fileinfo

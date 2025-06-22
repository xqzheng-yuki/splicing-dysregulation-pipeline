#path retrieval and directory logic

get_data_dir <- function(run_number) {
  search_directory <- "/mnt/gtklab01/xiaoqing"
  matching_folders <- list.files(path = search_directory, pattern = run_number, full.names = TRUE)
  data_dir <- glue::glue("{matching_folders}/star/result")
  info(logger,glue("You have set your data directory as {data_dir}"))
  return(data_dir)
}
get_bw_path <- function(path) {
  data_dir <- path
  bw_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                    group=c(control_group,treatment_group)),
                             tag=condition_tag) |>
    mutate(bw_path=glue::glue("{data_dir}/unmapped_CTX_{group}_{tag}.bw"))
  return(bw_fileinfo)
}
get_bam_path <- function(path) {
  data_dir <- path
  bam_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                     group=c(control_group,treatment_group)),
                              tag=condition_tag) |>
    mutate(bam_path=glue::glue("{data_dir}/unmappedAligned.sortedByCoord.out_CTX_{group}_{tag}.bam"))
  return(bam_fileinfo)
}
run_decoy_match <- list(
  "Run2" = "decoy3",
  "Run3" = "decoy5",
  "Run4" = "decoy6"
)

get_decoy_dir <- function(run_number) {
  search_directory <- "/mnt/gtklab01/xiaoqing/Z-decoy"
  decoy_number <- run_decoy_match[[run_number]]
  decoy_dir <- glue::glue("{search_directory}/{decoy_number}")
  info(logger,glue("You have set your decoy directory as {decoy_dir}"))
  return(decoy_dir)
}
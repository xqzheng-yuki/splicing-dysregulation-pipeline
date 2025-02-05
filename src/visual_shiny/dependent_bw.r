# File: dependent_bw.r
# Purpose: Clean up and wrap up all self-defined small function used in visualize_mapping.r
# Note: This code file is extracted from visualization_bw.qmd, please refer back for detailed information and explaination.

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
  search_directory <- "/mnt/gtklab01/xiaoqing/decoy"
  decoy_number <- run_decoy_match[[run_number]]
  decoy_dir <- glue::glue("{search_directory}/{decoy_number}")
  info(logger,glue("You have set your decoy directory as {decoy_dir}"))
  return(decoy_dir)
}

get_intron_bed <- function(path) {
  data_dir <- path
  intron_path <- paste0({data_dir},"/intronic.bed")
  if (!file.exists(intron_path)) {
    intron_path <- paste0({data_dir},"/intronic.pure.bed")
    }
  intron_data <- import.bed(intron_path)
  return(intron_data)
}
get_mashmap_bed <- function(path) {
  data_dir <- path
  if (str_sub(as.character(path), -1, -1)=="3") {
    mashmap_path <- paste0({data_dir},"/genome_found_sorted.bed")
  } else {
    mashmap_path <- paste0({data_dir},"/mashmap_intergenic.bed")
  }
  mashmap_data <- import.bed(mashmap_path)
  return(mashmap_data)
}

get_exon_bed <- function(path) {
  data_dir <- path
  exon_path <- paste0({data_dir},"/exon_sort.bed")
  if (!file.exists(exon_path)) {
    exon_path <- paste0({data_dir},"/exons.merged.bed")
  }
  exon_data <- import.bed(exon_path)
  return(exon_data)
}

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

import_bam <- function(file,bam_selection) {
  bam_param <- ScanBamParam(which = bam_selection)
  bam_data <- as(GenomicAlignments::readGAlignments(file, param = bam_param),"GRanges")
  bam_BRG <- BRGenomics::makeGRangesBRG(bam_data)
  keepStandardChromosomes(bam_param,pruning.mode = "coarse")
}

merge_treatment_bam <- function(selection_range,treatment_n,baminfo=bam_fileinfo) {
  baminfo <- dplyr::filter(baminfo,treatment==treatment_n)
  bam_data <- split(map(baminfo$bw_path, function(bam_path) {
    param <- ScanBamParam(which = selection_range)
    as(GenomicAlignments::readGAlignments(bam_path, param = param), "GRanges")
  }), baminfo$tag)
  bam_range <- map(bam_data, ~ {
                  mergeGRangesData(.x,exact_overlaps = TRUE)
    })
  return(bam_range)
}

generate_shades <- function(base_color) {
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
  # this if made for shiny app
  if (geneName == "Gene name") {
    gene_id_match <- ""
  } else {
  matches <- which(tolower(elementMetadata(geneRanges_GRCm39)$gene_name) == tolower(geneName))
  appr_matches <- agrep(geneName, 
                   elementMetadata(geneRanges_GRCm39)$gene_name, 
                   ignore.case = TRUE,
                   max.distance = 0)
  if (length(matches) == 0) {
    warn(logger,"There is no exact match found.")
    if (length(appr_matches) > 0) {
    fatal(logger,"There are more than one gene matched with your provided pattern.")
    gene_name_match <- paste(geneRanges_GRCm39[appr_matches]$gene_name,collapse = ", ")
    warn(logger,paste0(gene_name_match," have been found."))
    fatal(logger,"Please double check the gene name or provide a gene id starting with ENSMUSG.")
    return(invisible())
  }}
  gene_id_match <- geneRanges_GRCm39[matches]$gene_id
  info(logger,paste0("The gene id is ",gene_id_match))}
  return(gene_id_match)
}

is_gene_name <- function(gene_input) {
  # Assuming gene names are non-numeric and IDs are numeric
  # Adjust the logic based on your gene naming convention
  return(!grepl("^ENSMUSG\\d+\\.\\d+", gene_input))
}

get_ylim <- function(first_set,second_set,condition_tag = NULL) {
  ylim <- c(0,max(score(unlist(GRangesList(c(first_set,second_set))))))
  if (!is.null(condition_tag)) {
    ylim <- c(0,max(score(unlist(GRangesList(c(first_set[[condition_tag]],second_set[[condition_tag]]))))))
  }
  return(ylim)
}

enlarge_gr <- function(goi) {
  gr <- geneRanges_GRCm39[goi]
  embiggen_factor <- round(width(gr)*0.1)
  start(gr) <- start(gr) - embiggen_factor
  end(gr) <- end(gr) + embiggen_factor
  return(gr)
}


#' add_new_plots
#' update the state of the app with the most recent plot
#' @param state
#' @param classic_plot_func 
#' @param dataset_plot_func 
#'
#' @return state
#' @export
#'
#' @examples
add_new_plots <- function(state, classic_plot_func, dataset_plot_func) {
  state$results <- append(state$results, list(
    list(classic = classic_plot_func, dataset = dataset_plot_func)
  ))
  state$current_index <- length(state$results) # Set to the latest added
  # state$results[[state$current_index + 1]] <- recordPlot()
  info(logger,paste0("The state is ",length(state$results)))
  state
}
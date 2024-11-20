# File: dependent_bw.r
# Purpose: Clean up and wrap up all self-defined small function used in visualize_mapping.r
# Note: This code file is extracted from visualization_bw.qmd, please refer back for detailed information and explaination.


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
  info(logger,paste0("The gene id is ",gene_id_match))
  return(gene_id_match)
}

get_ylim <- function(first_set,second_set,condition_tag = NULL) {
  ylim <- c(0,max(score(unlist(GRangesList(c(first_set,second_set))))))
  if (!is.null(condition_tag)) {
    ylim <- c(0,max(score(unlist(GRangesList(c(first_set[[condition_tag]],second_set[[condition_tag]]))))))
  }
  return(ylim)
}
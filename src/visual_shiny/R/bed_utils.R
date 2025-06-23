# loading BED files and merging logic   

get_intron_bed <- function(path) {
  data_dir <- path
  intron_path <- paste0({data_dir},"/intronic.bed")
  if (!file.exists(intron_path)) {
    intron_path <- paste0({data_dir},"/intronic.pure.bed")
    }
  intron_data <- import.bed(intron_path)
  return(intron_data)
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

## Bigwig file
import_bigwig <- function(file,bw_selection) {
  bw_data <- BRGenomics::makeGRangesBRG(import.bw(file,selection=rtracklayer::BigWigSelection(bw_selection)))
  keepStandardChromosomes(bw_data,pruning.mode = "coarse")
}
merge_treatment_bw <- function(selection_range,treatment_n,bwinfo=bw_fileinfo) {
  bwinfo <- dplyr::filter(bwinfo,treatment==treatment_n)
  bw_data <- split(map(bwinfo$bw_path,import_bigwig,bw_selection=selection_range),
        bwinfo$tag)
  bw_range <- map(bw_data, ~ {
                  BRGenomics::mergeGRangesData(.x,exact_overlaps = TRUE)
    })
  return(bw_range)
}

## BAM file
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
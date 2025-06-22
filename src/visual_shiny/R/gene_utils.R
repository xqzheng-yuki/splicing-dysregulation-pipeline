# gene ID and name resolution
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
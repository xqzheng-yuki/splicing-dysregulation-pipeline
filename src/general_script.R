library(log4r)
library(tidyverse)
library(ggplot2)
tsv_path="/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/count"
treatment_group <- c(104,108,128,154)
control_group <- c(120,125,147,148)
tags=c("d","m1","m2")
get_tsv_path <- function(path) {
  data_dir <- path
  tsv_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                     group=c(control_group,treatment_group)),
                              tag=tags) |>
    mutate(tsv_path=glue::glue("{data_dir}/count_CTX_{group}_{tag}.tsv"))
  return(tsv_fileinfo)
}
get_gene_name <- function(geneID){
  matches <- which(elementMetadata(geneRanges_GRCm39)$gene_id == geneID)
  if (length(matches) == 0) {
    warn(logger,"There is no match found.")}
  gene_name_match <- geneRanges_GRCm39[matches]$gene_name
  # info(logger,paste0("The gene ",geneID," match is ",gene_name_match))
  return(gene_name_match)
}

tsv_info <- get_tsv_path(tsv_path)

tsv <- split(map(tsv_info$tsv_path, ~ {
  read_tsv(.x, col_names = c("gene","count"), show_col_types = FALSE)}),interaction(tsv_info$group,tsv_info$tag,sep="_"))

## clean gene column
for (sample in c(treatment_group,control_group)) {
  whole <- paste0(sample,"_d")
  # clean gene id suffix
  tsv[[whole]][[1]] <- subset(tsv[[whole]][[1]][str_detect(tsv[[whole]][[1]][["gene"]],"ENSMUSG"),])
  tsv[[whole]][[1]][["gene"]] <- sub("_intronic","",tsv[[whole]][[1]][["gene"]])
  # add gene name column
  tsv[[whole]][[1]] <- tsv[[whole]][[1]] %>% mutate(gene_name = map_chr(tsv[[whole]][[1]][["gene"]],get_gene_name))
  print(head(tsv[[whole]][[1]],1))
}

## add column for gene name
str_detect(tsv[1][[1]]$gene,"_intronic")
tsv[["104_d"]][[1]][["gene"]][1]
map_chr(tsv[["104_d"]][[1]][["gene"]],get_gene_name)
my_tibble <- my_tibble %>% mutate(new_col = some_value)
ggplot(data=de, aes(x=log2FoldChange, y=pvalue)) + geom_point()
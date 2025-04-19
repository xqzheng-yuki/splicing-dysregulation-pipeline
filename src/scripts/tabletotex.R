library(gt)
sig_gene_table <- significant_genes %>% 
  rowwise() |> 
  dplyr::mutate(gene_name=get_gene_name(chr)) %>% 
  dplyr::select(gene_name,log2FoldChange,pvalue,padj,baseMean) |>
  ungroup() |>
  dplyr::select(-pvalue) |>
  gt() |>
  fmt_number(columns=c(log2FoldChange,baseMean)) |>
  fmt_scientific(columns=c(padj)) |>
  as_latex()
cat(sig_gene_table)



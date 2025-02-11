library(tidyverse)
library(here)
library(DESeq2)
data_dir <- "/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/count"

treatment_group <- c(104,108,128,154)
control_group <- c(120,125,147,148)
tags=c("d","m1","m2")

tsv_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                     group=c(control_group,treatment_group)),
                              tag=tags) |>
    mutate(tsv_path=glue::glue("{data_dir}/count_CTX_{group}_{tag}.tsv"))
  
d_files <- 
intron_d_counts <- purrr::map(tsv_fileinfo |> filter(tag=='d') |> pull(tsv_path),read_tsv,col_names=c('chr','count')) |>
    purrr::reduce(inner_join,by='chr') |> 
    filter(str_starts(chr,"ENSMUSG"))

colnames(intron_d_counts)[-1] <- paste("CTX",tsv_fileinfo |> filter(tag=='d') |> pull(group),sep="_")

dds <- DESeqDataSetFromMatrix(column_to_rownames(intron_d_counts,'chr'),
    colData = tsv_fileinfo |> filter(tag=="d"),
    design = ~ treatment)

dds <- DESeq(dds)
res <- results(dds)
res_df <- lfcShrink(dds,coef=2,res=res)
data.frame(res_df) |> rownames_to_column('chr') |> tibble() |> arrange(padj)

## https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
ggplot(data=res_df, aes(x=padj, y=-log10(pvalue))) + geom_point()
p <- ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
res_df$diffexpressed <- "NO"
res_df$diffexpressed[res_df$log2FoldChange > 0.6 & res_df$pvalue < 0.05] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -0.6 & res_df$pvalue < 0.05] <- "DOWN"
p <- ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
res_df$delabel <- NA
res_df$delabel[res_df$diffexpressed != "NO"] <- rownames(res_df)[res_df$diffexpressed != "NO"]

ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

library(ggrepel)
ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

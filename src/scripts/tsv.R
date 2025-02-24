library(tidyverse)
library(here)
library(DESeq2)
library(ggrepel)
library(ggplot2)

source("/mnt/cbis/home/e1124735/Capstone/src/visual_shiny/dependent_2.R")
data_dir <- "/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/count"

tags=c("d","m1","m2")

# filter parameters
max_padj <- 0.05
min_lfc <- 0.5

# import data into dataframe
tsv_fileinfo <- expand_grid(tibble(treatment=rep(c("control","treatment"),each=4),
                                   group=c(control_group,treatment_group)),
                            tag=tags) |>
  mutate(tsv_path=glue::glue("{data_dir}/count_CTX_{group}_{tag}.tsv"))

intron_d_counts <- purrr::map(tsv_fileinfo |> dplyr::filter(tag=='d') |> pull(tsv_path),read_tsv,col_names=c('chr','count')) |>
  purrr::reduce(inner_join,by='chr') |> 
  dplyr::filter(str_starts(chr,"ENSMUSG"))

colnames(intron_d_counts)[-1] <- paste("CTX",tsv_fileinfo |> dplyr::filter(tag=='d') |> pull(group),sep="_")

# create deseq2 experiment
dds <- DESeqDataSetFromMatrix(column_to_rownames(intron_d_counts,'chr'),
                              colData = tsv_fileinfo |> dplyr::filter(tag=="d"),
                              design = ~ treatment)

dds <- DESeq(dds)
res <- results(dds)
res_df <- lfcShrink(dds,coef=2,res=res)
res_db <- data.frame(res_df) |> rownames_to_column('chr') |> tibble()

# subset of top 20 genes
significant_genes <- res_db %>%
  dplyr::filter(!is.na(log2FoldChange) & !is.na(padj) & padj > 0) %>%
  dplyr::filter(!is.na(padj) & padj < max_padj & (log2FoldChange > min_lfc | log2FoldChange < -min_lfc) ) %>%
  arrange(padj) %>% 
  slice_head(n = 20)
significant_genes$chr <- sub("_intronic","",significant_genes$chr)

## ref: https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

# add differentiated expressed and label col
res_db$diffexpressed <- NA
res_db$diffexpressed[res_db$log2FoldChange > min_lfc & res_db$padj < max_padj] <- "UP"
res_db$diffexpressed[res_db$log2FoldChange < -min_lfc & res_db$padj < max_padj] <- "DOWN"
table(res_db$diffexpressed,useNA = "ifany")

res_db$delabel <- NA
res_db <- res_db |> 
  rowwise() |>
  mutate(chr=str_extract(chr,"ENSMUSG\\d+.\\d+"),
         delabel=ifelse(padj < max_padj & abs(log2FoldChange) >= min_lfc,map_chr(chr,get_gene_name),NA))

# plot `p`: decoy data in volcano plot version
p <- ggplot() +
  geom_point(
    data=res_db,
    mapping = aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed),
    stat = 'identity',) + 
  theme_minimal() +
  geom_label_repel(data = res_db,mapping = aes(x=log2FoldChange, y=-log10(padj),label=delabel)) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("blue", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed") +
  ggtitle("Decoy Result as base")

# import spliceat golden standard for comparison
list_ctx<-"/mnt/gtklab01/xiaoqing/golden_standard/tdp43_nestin_ctx_e14_gene_mode_lrt_wald_join.csv"
golden_df <- read.csv(list_ctx, sep=",")

g_res <- tibble(
  log2FoldChange = golden_df$b_wald,
  pvalue = golden_df$pval_wald,
  padj = golden_df$qval_wald,
  chr = golden_df$target_id) |>
  mutate(delabel=ifelse(padj < max_padj & abs(log2FoldChange) > min_lfc, map_chr(chr,get_gene_name),NA)) |>
  relocate(chr)

# add differentiated expressed feature for later color diff
g_res$diffexpressed <- NA
g_res$diffexpressed[g_res$log2FoldChange > min_lfc & g_res$padj < max_padj] <- "UP"
g_res$diffexpressed[g_res$log2FoldChange < -min_lfc & g_res$padj < max_padj] <- "DOWN"

# plot `p2`: spliceat-filtered gene list in volcano plot version
p2 <- ggplot() +
  geom_point(
    data=g_res,
    mapping = aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed),
    stat = 'identity',) + 
  theme_minimal() +
  geom_label_repel(data = g_res,mapping = aes(x=log2FoldChange, y=-log10(padj),label=delabel)) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("blue", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed") +
  ggtitle("Spliceat gold standard as base")

# overlap two datasets for overlay analysis
overlap <- inner_join(res_db,g_res,by='delabel',na_matches="never") |>
  dplyr::select(chr=chr.x,FC_decoy=log2FoldChange.x,PV_decoy=padj.x,FC_spl=log2FoldChange.y,PV_spl=padj.y,delabel)
# additional layer of points for overlap indication
p3 <- geom_point(data = overlap, mapping = aes(x=FC_decoy,y=-log10(PV_decoy)), shape=11, size=2) # for decoy set
p4 <- geom_point(data = overlap, mapping = aes(x=FC_spl,y=-log10(PV_spl)), shape=11, size=2) # for spliceat set

# Whole graph of decoy set of differential gene with overlap gene shown as stared point
p+p3
# Whole graph of spliceat augmented cryptic gene with overlap gene shown as stared point
p2+p4
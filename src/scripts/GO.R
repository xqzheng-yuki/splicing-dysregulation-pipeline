library(clusterProfiler)
if (!require("org.Mm.eg.db", quietly = TRUE))
  install.packages("org.Mm.eg.db")
library(org.Mm.eg.db)

# note the variable is from tsv.R

gene_entrez <- bitr(all_significant$chr,
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
go_enrich <- enrichGO(gene = gene_entrez$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
dotplot(go_enrich)
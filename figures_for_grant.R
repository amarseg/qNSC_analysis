library(tidyverse)
library(clusterProfiler)
source('useful_funs.R')
library(DESeq2)
library(org.Mm.eg.db)

pval_thr <- 0.05
fold_change_thr <- 1.5

coculture_data <- read_delim('../supp1.txt', delim = '\t') %>%
  select(Gene_names, log2FoldChange_coculture_vs_control:padj_coculture_vs_control) %>%
  filter(padj_coculture_vs_control < pval_thr & abs(log2FoldChange_coculture_vs_control) > fold_change_thr) %>%
  add_column(direction = ifelse(.$log2FoldChange_coculture_vs_control > 0,'up','down'))

entrez_ids <- bitr( geneID = coculture_data$Gene_names, 
                                  fromType = 'ALIAS', 
                                  toType = 'ENTREZID',
                                  OrgDb = org.Mm.eg.db)

coculture_data <- left_join(coculture_data, entrez_ids, by = c('Gene_names' = 'ALIAS'))

go <- compareCluster(coculture_data, ENTREZID ~ direction, fun = 'enrichGO', OrgDb = org.Mm.eg.db)
kegg <- compareCluster(coculture_data, ENTREZID ~ direction, fun = 'enrichKEGG', org = 'mmu')
david <- compareCluster(coculture_data, ENTREZID ~ direction, fun = 'enrichDAVID',
                        david.user = 'a.martinez-segura13@imperial.ac.uk')

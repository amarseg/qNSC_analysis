source('useful_funs.R')
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
source('useful_funs.R')
library(org.Mm.eg.db)
library(ggrepel)
library('wesanderson')

pval_thr <- 0.05
fold_change_thr <- 1

coculture_data <- read_delim('../supp1.txt', delim = '\t') %>%
  dplyr::select(Gene_names, log2FoldChange_coculture_vs_control:padj_coculture_vs_control) %>%
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


citokines <- read_delim('../IPA_files/cytokines.txt', delim = '\t', skip = 1)
ids <- strsplit(citokines$`Target molecules in dataset`, split = ',', fixed = T)
citokines$overlap <- unlist(lapply(ids, length))
citokines$mean_xp <- unlist(lapply(ids, summarise_exp,expr = as.data.frame(coculture_data[,c(1,2)])))

ggplot(citokines, 
       aes(x = mean_xp, 
           y = -log10(`p-value of overlap`), 
           label = `Upstream Regulator`, 
           colour = `Predicted Activation State`)) +
  geom_point(data = filter(citokines, !is.na(`Predicted Activation State`))) + 
  geom_text_repel(size = 2.5) + 
  theme_light() + 
  scale_colour_manual(values = wes_palette('BottleRocket2')[c(2,3)], name = 'Activation state')

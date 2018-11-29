source('useful_funs.R')
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
source('useful_funs.R')
library(org.Mm.eg.db)
library(ggrepel)
library('wesanderson')
library('pheatmap')
library(ggpubr)

pval_thr <- 0.05
fold_change_thr <- 1

#####################Read data############

coculture_data_de <- read_delim('../supp1.txt', delim = '\t') %>%
  dplyr::select(Gene_names, log2FoldChange_coculture_vs_control:padj_coculture_vs_control, `Codega et al.`,`Boadilla et al.`) %>%
  filter(padj_coculture_vs_control < pval_thr & abs(log2FoldChange_coculture_vs_control) > fold_change_thr) %>%
  add_column(direction = ifelse(.$log2FoldChange_coculture_vs_control > 0,'up','down'))

entrez_ids <- bitr( geneID = coculture_data_de$Gene_names, 
                                  fromType = 'ALIAS', 
                                  toType = 'ENTREZID',
                                  OrgDb = org.Mm.eg.db)

##################Enrichments using clusterProfiler#####################
coculture_data_de <- left_join(coculture_data_de, entrez_ids, by = c('Gene_names' = 'ALIAS'))

#go <- compareCluster(data = coculture_data, ENTREZID ~ direction, fun = 'enrichGO', OrgDb = org.Mm.eg.db)
# dotplot(go, showCategory = 20)
kegg <- compareCluster(data = coculture_data_de, ENTREZID ~ direction, fun = 'enrichKEGG', org = 'mmu')
dotplot(kegg, showCategory = 20)
david <- compareCluster(data = coculture_data_de, ENTREZID ~ direction, fun = 'enrichDAVID',david.user = 'a.martinez-segura13@imperial.ac.uk')
dotplot(david, showCategory = 20)

#########################Bubble plot with cytokines that IPA predicts#####################

citokines <- read_delim('../IPA_files/cytokines.txt', delim = '\t', skip = 1)
ids <- strsplit(citokines$`Target molecules in dataset`, split = ',', fixed = T)
citokines$overlap <- unlist(lapply(ids, length))
citokines$mean_xp <- unlist(lapply(ids, summarise_exp,expr = as.data.frame(coculture_data_de[,c(1,2)])))

ggplot(data = filter(citokines, !is.na(`Predicted Activation State`)), 
       aes(x = mean_xp, 
           y = -log10(`p-value of overlap`), 
           label = `Upstream Regulator`, 
           colour = `Predicted Activation State`,
           size = overlap)) +
  geom_point() + 
  geom_text_repel(size = 3) + 
  theme_light()


##############Barplot using the cytokines from MGI#################

mgi_cito <- read_delim('cytokines_from_MGI.txt', delim = '\t')
mgi_cito <- unique(mgi_cito$Symbol)

citos <- filter(coculture_data_de, Gene_names %in% mgi_cito)

ggplot(citos, aes(x = reorder(Gene_names,-log2FoldChange_coculture_vs_control), y = log2FoldChange_coculture_vs_control)) +
  geom_bar(stat = 'identity')

cocultures_all <- read_delim('../supp1.txt', delim = '\t') %>%
  dplyr::select(Gene_names, log2FoldChange_coculture_vs_control:padj_coculture_vs_control)

all_citos <- filter(cocultures_all, Gene_names %in% mgi_cito)

ggplot(all_citos, aes(x = reorder(Gene_names,-log2FoldChange_coculture_vs_control), y = log2FoldChange_coculture_vs_control)) +
  geom_bar(stat = 'identity')

###########################Look at the targets from IFNG######################
ifng_targets <- filter(citokines, `Upstream Regulator` == 'IFNG') %>%
  dplyr::select(`Target molecules in dataset`)

ifng_targets <- str_split(ifng_targets, pattern = ',')[[1]] %>%
  str_lowercase() %>%
  str_ucfirst() %>%
  as.data.frame()

ifng_entrez <- bitr(geneID = ifng_targets[,1], fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
ifng_targets <- left_join(ifng_targets, ifng_entrez, by = c('.'='ALIAS'))

ifng_exp <- inner_join(ifng_targets, coculture_data_de, by = c('.'='Gene_names'))

s <- dplyr::select(ifng_exp,ENTREZID.x, direction) %>%
    split(ifng_exp, f = ifng_exp$direction, drop = T)

final <- list(s[[1]]$ENTREZID.x, s[[2]]$ENTREZID.x)
names(final) <- c('down','up')

go_2 <- compareCluster(final, fun = 'enrichGO', OrgDb = org.Mm.eg.db)
kegg_2 <- compareCluster(final, fun = 'enrichKEGG', org = 'mmu')
david_2 <- compareCluster(final, fun = 'enrichDAVID',
                        david.user = 'a.martinez-segura13@imperial.ac.uk')
dotplot(kegg_2, showCategory = 20)

########################How many of these guys are in Codega and Bobadilla#######################

ggplot(filter(citos, !is.na(`Codega et al.`)), aes(x = `Codega et al.`)) +
  geom_bar()

ggplot(filter(citos, !is.na(`Boadilla et al.`)), aes(x = `Boadilla et al.`)) +
  geom_bar()

ggplot(citos, aes(x = interaction(`Boadilla et al.`, `Codega et al.`))) +
  geom_bar() 
#############################IFNGN targets according to IPA############################
ipa_targets <- read_csv('ifng_targets_IPA.txt') %>%
  colnames() %>%
  str_lowercase() %>%
  str_ucfirst() %>%
  as.data.frame() 

ipa_targets_df <- add_column(cocultures_all, target = ifelse(cocultures_all$Gene_names %in% ipa_targets[,1]
                                                             ,'yes','no'))

ggplot(ipa_targets_df, aes( y = log2FoldChange_coculture_vs_control, x = target)) +
  geom_boxplot(notch = T) +
  stat_compare_means(method = "wilcox.test")

#####################codega citokines################################################
citokines_codega <- read_delim('cytokines_codega.txt', delim = '\t', skip = 1)
ids <- strsplit(citokines_codega$`Target molecules in dataset`, split = ',', fixed = T)
citokines_codega$overlap <- unlist(lapply(ids, length))
citokines_codega$mean_xp <- unlist(lapply(ids, summarise_exp,expr = as.data.frame(coculture_data_de[,c(1,2)])))

all_citokines <- inner_join(citokines, citokines_codega, by = 'Upstream Regulator')

ggplot(data = all_citokines, 
       aes(x = mean_xp.x, 
           y = -log10(`p-value of overlap.x`), 
           label = `Upstream Regulator`, 
           colour = `Predicted Activation State.y`,
           size = overlap.x)) +
  geom_point() + 
  facet_wrap(~`Predicted Activation State.x`) +
  geom_text_repel(size = 3) + 
  theme_light()

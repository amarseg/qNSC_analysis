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

entrez_ids <- bitr(geneID = coculture_data_de$Gene_names, 
                                  fromType = 'ALIAS', 
                                  toType = 'ENTREZID',
                                  OrgDb = org.Mm.eg.db)

##################Enrichments using clusterProfiler#####################
coculture_data_de <- inner_join(coculture_data_de, entrez_ids, by = c('Gene_names' = 'ALIAS'))

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

#go_2 <- compareCluster(final, fun = 'enrichGO', OrgDb = org.Mm.eg.db)
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

###############################Volcano plots of cytokines and chemokines########################
chemokines <- read_delim('chemokine_Activity_mgi.txt', delim = '\t')
chemokines <- unique(chemokines$Symbol)

chemokines = data.frame(gene_name = chemokines, type = 'Chemokines')
cytokines = data.frame(gene_name = mgi_cito, type = 'Cytokines')

all_types <- bind_rows(chemokines, cytokines)

coculture_annotated <- cocultures_all %>%
  left_join(all_types, by = c('Gene_names' = 'gene_name'))

kines_only <- filter(coculture_annotated, !is.na(type))

ggplot(kines_only, aes(x = log2FoldChange_coculture_vs_control, 
                                y = -log10(padj_coculture_vs_control),
                                colour = type,
                                label = Gene_names)) +
  geom_point(data = subset(kines_only, -log10(padj_coculture_vs_control) < 2)) +
  geom_text(data = subset(kines_only, -log10(padj_coculture_vs_control) > 2)) +
  facet_wrap(~type) +
  theme_light() +
  scale_colour_manual(values = wes_palette('Moonrise2'), name = 'Molecule type') +
  geom_hline(yintercept=2, linetype = 'dashed', colour = 'darkgray')

#######################GSEA#################################################
entrez_ids <- bitr( geneID = cocultures_all$Gene_names, 
                    fromType = 'ALIAS', 
                    toType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db)

##################Enrichments using clusterProfiler#####################


pathways <- gmtPathways('h.all.v6.2.symbols.gmt')


## feature 1: numeric vector
geneList = cocultures_all$log2FoldChange_coculture_vs_control
## feature 2: named vector
names(geneList) = as.character(cocultures_all$Gene_names)
names(geneList) <- toupper(names(geneList))
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

fgseaRes <- fgsea(pathways = pathways, 
                  stats = geneList,
                  minSize=10,
                  maxSize=1000,
                  nperm=100000)


save(fgseaRes, file  =  'gsea_coculture_results.Rdata')

sig_Path=fgseaRes[which(fgseaRes[,"padj"] < 0.05),]
sig_Path=sig_Path[order(sig_Path[,"NES"],decreasing=T),]

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15),]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], geneList, fgseaRes, gseaParam = 0.5)

plotEnrichment(pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
               geneList) + labs(title="IFNG response")

plotEnrichment(pathways[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
               geneList) + labs(title="IFNA response")
up <- fgseaRes[ES > 0]

ggplot(top_n(fgseaRes, n = 15, padj), aes( x = -log10(padj), y = size, label = pathway)) +
  geom_point() +
  geom_text_repel()

#######################################Merge codega with coculture dataset##################
aNSC_codega <- read_csv('../Codega_supp/aNSC_signature.csv', skip = 2) %>%
  dplyr::select(`Gene Symbol`, `Corrected p-value`, `Fold Change`) %>%
  add_column(cells = 'aNSC')
qNSC_codega <- read_csv('../Codega_supp/qNSC_signature.csv', skip = 2) %>%
  dplyr::select(`Gene Symbol`, `Corrected p-value`, `Fold Change`) %>%
  add_column(cells = 'qNSC')

codega <- bind_rows(aNSC_codega, qNSC_codega)

coculture_codega <- left_join(cocultures_all, codega, by = c('Gene_names' = 'Gene Symbol')) %>%
  dplyr::select(-ENTREZID.x, -ENTREZID.y) %>%
  write_csv('coculture_data_plus_codega.csv')

coculture_codega_annot <- inner_join(coculture_codega, all_types , by = c('Gene_names' = 'gene_name'))

pheatmap(dplyr::select(coculture_codega_annot, log2FoldChange_coculture_vs_control, 'Fold Change'))

true_ones <- inner_join(codega, all_types, by = c('Gene Symbol' = 'gene_name'))
filter(coculture_codega_annot, type == 'Chemokines')
filter(coculture_codega_annot, type == 'Cytokines')

######################################DO GSEA with codega dataset###############################
codega <- read_csv('../Codega_supp/all_codega.csv')
## feature 1: numeric vector
geneList = codega$`Fold Change`
## feature 2: named vector
names(geneList) = as.character(codega$`Gene Symbol`)
names(geneList) <- toupper(names(geneList))
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

fgseaRes <- fgsea(pathways = pathways, 
                  stats = geneList,
                  minSize=10,
                  maxSize=1000,
                  nperm=100000)


save(fgseaRes, file  =  'gsea_codega_results.Rdata')

sig_Path=fgseaRes[which(fgseaRes[,"padj"] < 0.05),]
sig_Path=sig_Path[order(sig_Path[,"NES"],decreasing=T),]

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15),]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[[topPathways]], geneList, fgseaRes, gseaParam = 0.5)

plotEnrichment(pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
               geneList) + labs(title="IFNG response")

plotEnrichment(pathways[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
               geneList) + labs(title="IFNA response")
up <- fgseaRes[ES > 0]

ggplot(top_n(fgseaRes, n = 15, padj), aes( x = -log10(padj), y = size, label = pathway)) +
  geom_point() +
  geom_text_repel()
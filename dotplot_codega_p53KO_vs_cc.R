##############figure 5 according to Tim's draft
rm(list = ls())
########## Libraries and setwd#####
#Libraries

library('DESeq2')
library('pheatmap')
library('IHW')
library('clusterProfiler')
library('org.Mm.eg.db')
library('RDAVIDWebService')
library('RColorBrewer')
library('tidyverse')
library('ggrepel')
setwd('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/figure_2')
source('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/Grouping.R')
source('useful_funs.R')
library('gridExtra')

#####
counts <- read.delim('../supp1.txt', header = T, strings = F, na.strings = '<NA>')

dataset <- counts[,c(1,5,6,7,20,21)]
dataset$plot <- -log10(as.numeric(dataset$padj_p53KO_coculture_vs_flox_coculture))
############Volcano plot with colored codega

special.boadilla <- dataset$Boadilla.et.al. != 'NA'
p <- ggplot(dataset, aes(x = log2FoldChange_p53KO_coculture_vs_flox_coculture, y = plot))
p2<-p + geom_point(alpha = 0.20, color = 'darkgrey') + theme_light() + geom_point(data=dataset[special.boadilla,], aes(fill=Boadilla.et.al.), color="black", alpha=0.8, size=2, shape=21)

special.codega <- dataset$Codega.et.al. != 'NA'
p <- ggplot(dataset, aes(x = log2FoldChange_p53KO_coculture_vs_flox_coculture, y = plot))
p1<- p + geom_point(alpha = 0.20, color = 'darkgrey')+ theme_light() + geom_point(data=dataset[special.codega,], aes(fill=Codega.et.al.), color="black", alpha=0.8, size=2, shape=21)

grid.arrange(p1,p2, ncol = 2)

#################who is up and belongs to codega and boadilla
sig_up <- dataset[dataset$padj_p53KO_coculture_vs_flox_coculture < 0.05 & dataset$log2FoldChange_p53KO_coculture_vs_flox_coculture > 0,]
only_physio_up <- sig_up[sig_up$Codega.et.al. == 'aNSC' & sig_up$Boadilla.et.al. == 'aNSC',]
sig_up_ent <- bitr(only_physio_up$Gene_names, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

kegg_up_aNSC <- enrichDAVID(sig_up_ent$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk')
dotplot(kegg_up_aNSC)

kegg_up_aNSC <- enrichKEGG(sig_up_ent$ENTREZID, organism = 'mmu')
dotplot(kegg_up_aNSC)

#################who is down and belongs to codega and boadilla
sig_down <- dataset[dataset$padj_p53KO_coculture_vs_flox_coculture < 0.05 & dataset$log2FoldChange_p53KO_coculture_vs_flox_coculture < 0,]
only_physio_down <- sig_down[sig_down$Codega.et.al. == 'qNSC' | sig_down$Boadilla.et.al. == 'qNSC',]
sig_down_ent <- bitr(only_physio_down$Gene_names, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

kegg_down_aNSC <- enrichDAVID(sig_down_ent$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk')
dotplot(kegg_down_aNSC)
kegg_df <- as.data.frame(kegg_down_aNSC)

down_kegg <- enrichKEGG(sig_down_ent$ENTREZID, organism = 'mmu')
dotplot(down_kegg)
ppar <- down_kegg[grep(down_kegg$Description, pattern = 'PPAR', fixed = T),]
ppar_names <- as.data.frame(strsplit(ppar$geneID, split = '/', fixed = T))
ppar_names <- bitr(ppar_names[,1], fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = org.Mm.eg.db)

fatty <- kegg_df[grep(kegg_df$Description, pattern = 'fatty', fixed = T),]
fatty <- as.data.frame(strsplit(fatty$geneID, split = '/', fixed = T))
t <- bitr(fatty[,1], fromType = 'ENTREZID', toType = 'ALIAS', OrgDb = org.Mm.eg.db)
write.table(t,'physio_down_regulated.csv', sep = ',')

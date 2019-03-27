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

gmtfile <- system.file("extdata", "c7.all.v6.2.entrez.gmt", package="clusterProfiler")
c7 <- read.gmt(gmtfile)

coculture_data_de <- read_delim('../supp1.txt', delim = '\t') %>%
  dplyr::select(Gene_names, log2FoldChange_coculture_vs_control:padj_coculture_vs_control, `Codega et al.`,`Boadilla et al.`) %>%
  filter(padj_coculture_vs_control < pval_thr & abs(log2FoldChange_coculture_vs_control) > fold_change_thr) %>%
  add_column(direction = ifelse(.$log2FoldChange_coculture_vs_control > 0,'up','down'))

library(fgsea)
#
ranked_CC_BB=DATA[which(DATA[,7] < 0.05),1]
names(ranked_CC_BB)=row.names(DATA)[which(DATA[,7] < 0.05)]
ranked_CC_BB=sort(ranked_CC_BB,decreasing=T)
#
ranked_CC_BB=DATA[,1]
names(ranked_CC_BB)=row.names(DATA)
ranked_CC_BB=sort(ranked_CC_BB,decreasing=T)
#
#a1=ranked_CC_BB[which(ranked_CC_BB >= 1)]
#a2=ranked_CC_BB[which(ranked_CC_BB <= -1)]
#slim_CC_BB=c(a1,a2)
#slim_CC_BB=sort(slim_CC_BB,decreasing=T)


fgseaRes <- fgsea(pathways = Purged_lists[test], 
                  stats = ranked_CC_BB,
                  minSize=10,
                  maxSize=1000,
                  nperm=100000)

head(fgseaRes[order(pval),1:5],n=10)
tail(fgseaRes[order(pval),1:5],n=20)
sig_Path=fgseaRes[which(fgseaRes[,"padj"] < 0.05),]
sig_Path=sig_Path[order(sig_Path[,"NES"],decreasing=T),]
#
i=1
plotEnrichment(Purged_lists[fgseaRes[order(pval)[i],pathway]],
               ranked_CC_BB) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)
#
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(Purged_lists[topPathways], ranked_CC_BB, fgseaRes, gseaParam = 0.5)
plotGseaTable(Purged_lists, slim_CC_BB, fgseaRes, gseaParam = 0.5)
plotGseaTable(Purged_lists[as.character(sig_Path[,pathway])], ranked_CC_BB, sig_Path, gseaParam = 0.5)
fgseaRes1=fgseaRes[order(fgseaRes[,"NES"],decreasing = T),]
plotGseaTable(Purged_lists[as.character(unlist(fgseaRes1[,1]))], ranked_CC_BB, fgseaRes1, gseaParam = 0.5)

##########
#I take the upstream regulator table for upregulated, downregulated in coculture
#Same analysis was performed with aNSC and qNSC lists
#IPA analysis carried out taking into account the fold change!!! not only the overlap
#Take only TFs in common!!!!

library('tidyverse')
library('pheatmap')
library('wesanderson')
library('RColorBrewer')
source('../../figure_2/useful_funs.R')
setwd('P:/Simona_project/IPA_files/TF_files/')

aNSC <- read.delim('aNSC_fold_change.txt', header = T, strings = F, skip = 2)
aNSC <- aNSC[,c(1,7)]
aNSC <- aNSC[aNSC$p.value.of.overlap < 0.05,]
qNSC <- read.delim('qNSC_fold_change.txt', header = T, strings = F, skip = 2)
qNSC <- qNSC[,c(1,7)]
qNSC <- qNSC[qNSC$p.value.of.overlap < 0.05,]
up <- read.delim('up_regulated_tfs.txt', header = T, strings = F, skip = 2)
up <- up[,c(1,7)]
up <- up[up$p.value.of.overlap < 0.05,]
down <- read.delim('down_regulated_tfs.txt', header = T, strings = F, skip = 2)
down <- down[,c(1,6)]
down <- down[down$p.value.of.overlap < 0.05,]

down_aNSC <- inner_join(aNSC, down, by = 'Upstream.Regulator')
row.names(down_aNSC) <- down_aNSC$Upstream.Regulator
down_aNSC <- down_aNSC[,-1]
pval_trans <- -log10(down_aNSC)
colnames(pval_trans) <- c('Down_coculture','aNSC')

cutoff <- 5
pretty_col <- brewer.pal(9,'PuBuGn')
pval_trans <- pval_trans[which(pval_trans[,1] > cutoff | pval_trans[,2] > cutoff),]
pheatmap(pval_trans, cluster_cols = F, color = pretty_col, display_numbers = T, cex = 0.75)

up_qNSC <- inner_join(qNSC, up , by = 'Upstream.Regulator')
row.names(up_qNSC) <- up_qNSC$Upstream.Regulator
up_qNSC <- up_qNSC[,-1]
pval_trans <- -log10(up_qNSC)
colnames(pval_trans) <- c('up_coculture','qNSC')
pval_trans <- pval_trans[which(pval_trans[,1] > cutoff | pval_trans[,2] > cutoff),]

pheatmap(pval_trans, cluster_cols = F, color = pretty_col, display_numbers = T, cex = 0.75)

go_fatty_acid_degradation<- read.delim('../../go_fatty_Acid.txt', header = F, strings = F)
fatty_ids <- unique(go_fatty_acid_degradation[,2])
fatty_ids <- append(fatty_ids, 'Acad11')

heatmapper(fatty_ids)

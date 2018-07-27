rm(list = ls())
########## Libraries and setwd#####
#Libraries

library('DESeq2')
library('pheatmap')
library('clusterProfiler')
library('org.Mm.eg.db')
library('RDAVIDWebService')
library('RColorBrewer')
library('tidyverse')
library('pathview')
setwd('P:/Simona_project/figure_1')
source('../figure_2/useful_funs.R')
######Load data for coculture
supp <- read.delim('../supp1.txt', header = T, strings= F, na.strings = 'NA')
supp[supp=='<NA>'] <- NA

cc_data <- supp[,c(1:4, 20,21)]
entrez_ids <- bitr(cc_data$Gene_names, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

cc_data <- inner_join(cc_data, entrez_ids, by = c('Gene_names' = 'SYMBOL'))
cc_data <- cc_data[which(!duplicated(cc_data$Gene_names)),]

de_obj <- cc_data[,c(2:4)]
row.names(de_obj) <- cc_data$Gene_names

colnames(de_obj) <- c('log2FoldChange','stat','p.adjust' )
##############What is missing in vitro relative to in vivo?

cc_data$DE <- 'Not DE'
cc_data[which(cc_data$padj_coculture_vs_control < 0.05),]$DE <- 'DE'

missing_in_vitro <- cc_data[which(cc_data$DE == 'Not DE' & (!is.na(cc_data$Codega.et.al.) | !(is.na(cc_data$Boadilla.et.al.)) )),]

table(missing_in_vitro[,5:6], useNA = 'always') ##Some guys are not in common (actually many!)

##For now we only use codega
qNSC_missing <- missing_in_vitro[which(missing_in_vitro$Codega.et.al. == 'qNSC'),]
aNSC_missing <- missing_in_vitro[which(missing_in_vitro$Codega.et.al. == 'aNSC'),]

kegg_qNSC <- enrichKEGG(qNSC_missing$ENTREZID, organism = 'mmu'))
kegg_aNSC <- enrichKEGG(aNSC_missing$ENTREZID, organism = 'mmu')

qNSC_bub <- i_plot_bubbles(de_obj, enrichment_obj = kegg_qNSC@result, labels = 0, ID = T,title = 'test')
aNSC_bub <- i_plot_bubbles(de_obj, enrichment_obj = kegg_aNSC@result, labels = 0, ID = T,title = 'test')

all_bub <-join_bubble_objects(qNSC_bub, aNSC_bub)
all_bub$type <- 'aNSC'
all_bub[which(all_bub$term %in% kegg_qNSC@result$Description),]$type <- 'qNSC'

p <- ggplot(all_bub, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel(size = 4, color = 'black') + theme_light() + xlim(-2,2) 

res_qNSC <- kegg_qNSC@result
res_qNSC$Cell.type <- 'qNSC'
res_aNSC <- kegg_aNSC@result
res_aNSC$Cell.type <- 'aNSC'

all_res <- rbind(res_qNSC, res_aNSC)

split_ids <- strsplit(all_res$geneID, split = '/', fixed = T)
split_symbols <- lapply(split_ids, FUN = bitr, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = org.Mm.eg.db)
split_symb <- lapply(split_symbols, function(x){x$SYMBOL})
col_symb <- as.vector(lapply(split_symb, FUN = paste, collapse = '/'))
all_res$geneID <- unlist(col_symb)

test <- inner_join(all_res, all_bub[,1:2], by = c('Description' = 'term'))
colnames(test)[11] <- 'Mean log2FoldChange'
write.table(test, 'enrichment_missing_invivo.txt', sep = '\t')

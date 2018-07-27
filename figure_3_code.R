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
setwd('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/figure_3')
source('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/Grouping.R')
source('../figure_2/useful_funs.R')

#############Load DESeq object and extract BMP vs no BMP comparison
load('../figure_2/grouped_DESeq.RData')

bmp_de <- results(grouped_ds, contrast=c('group','floxbmp','floxcontrol_bmp'))

pval_threshold <- 0.05
fc_threshold <- 0

write.table(bmp_de, 'bmp_deseq.csv', sep =',')

bmp_sig <- bmp_de[bmp_de$padj < pval_threshold & !is.na(bmp_de$padj),]
bmp_sig <- bmp_sig[abs(bmp_sig$padj) > fc_threshold,]
bmp_sig <- as.data.frame(bmp_sig)

up_bmp <- bmp_sig[bmp_sig$log2FoldChange > 0,]
up_id <- bitr(row.names(up_bmp), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
down_bmp <- bmp_sig[bmp_sig$log2FoldChange < 0,]
down_id <- bitr(row.names(down_bmp), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

###########################Do all the enrichments!!!!########################
##KEGG
pdf('bmp_DE_enrichments.pdf', width = 25, height = 25)
up_kegg <- enrichKEGG(up_id$ENTREZID, organism = 'mmu')
enrichMap(up_kegg)
down_kegg <- enrichKEGG(down_id$ENTREZID, organism = 'mmu')
enrichMap(down_kegg) #p53 signalling pathway!

###GO
up_GO <- enrichGO(up_id$ENTREZID, OrgDb = org.Mm.eg.db,keytype = 'ENTREZID')
enrichMap(up_GO)
down_GO <- enrichGO(down_id$ENTREZID, OrgDb = org.Mm.eg.db,keytype = 'ENTREZID')
enrichMap(down_GO)

####DAVID
up_DAVID <- enrichDAVID(up_id$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk')
enrichMap(up_DAVID)
down_DAVID <- enrichDAVID(down_id$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk')
enrichMap(down_DAVID)
dev.off()
####################################Comparison with genes DE in coculture#####################
cc_de <- results(grouped_ds, contrast=c('group','floxendo','floxcontrol_endo'))

cc_sig <- cc_de[cc_de$padj < pval_threshold & !is.na(cc_de$padj),]
cc_sig <- cc_sig[abs(cc_sig$padj) > fc_threshold,]
cc_sig <- as.data.frame(cc_sig)

cc_bmp_common <- intersect(row.names(cc_sig), row.names(bmp_sig))
common_ids <- bitr(cc_bmp_common, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db) 
common_kegg <- enrichKEGG(common_ids$ENTREZID, org = 'mmu')
enrichMap(common_kegg)

#############################Heatmap of genes in common#########################
cc_common <- cc_sig[which(row.names(cc_sig) %in% common_ids$ALIAS),]
bmp_common <- bmp_sig[which(row.names(bmp_sig) %in% common_ids$ALIAS),]
pdf('bmp_DE_in_cc.pdf')
all_common <- merge(cc_common, bmp_common, by = 'row.names')
plot(x = all_common$log2FoldChange.x, y = all_common$log2FoldChange.y)
#Some of the response is in the opposite direction?
#who are these guys

up_cc_down_bmp <- all_common[which(all_common$log2FoldChange.x > 0 & all_common$log2FoldChange.y < 0 ),]
up_down_id <- bitr(up_cc_down_bmp$Row.names, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
down_cc_up_bmp <- all_common[which(all_common$log2FoldChange.x < 0 & all_common$log2FoldChange.y > 0 ),]
down_up_id <- bitr(down_cc_up_bmp$Row.names, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

good_list <- list(up_down = up_down_id$ENTREZID, down_up = down_up_id$ENTREZID)
test <- compareCluster(geneClusters = good_list, fun = 'enrichKEGG', org = 'mmu')
dotplot(test)

#################################nice bublle plot#########################
up_bub <- i_plot_bubbles(up_kegg, bmp_sig, 5, T, 'up') 
down_bub <- i_plot_bubbles(down_kegg, bmp_sig, 5, T, 'down') 

down_bub <- i_plot_bubbles(rbind(as.data.frame(down_kegg), as.data.frame(up_kegg)), bmp_sig, 5, T, 'down') 


all_bub <- join_bubble_objects(up_bub, down_bub)
p <- ggplot(all_bub, aes(x = means, y = -log10(pvalue), size = counter, fill = means, label = term))
p2<- p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel(size = 4, color = 'black') + theme_light() + scale_fill_distiller(palette = 'Spectral')
dev.off()
ggsave(p2, 'bubble_plot_DE_bmp.pdf')

comparison_and_enrichment(bmp_de, cc_de, labels = c('bmp','coculture'), plot_title = 'bmp_vs_coculture.pdf', k = 8)

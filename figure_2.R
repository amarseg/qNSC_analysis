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
setwd('P:/Simona_project/figure_2')
source('P:/Simona_project/06_Analysis/02_Data/Grouping.R')
source('useful_funs.R')
###############Load data######
#Load data
load('../DESeq_filtered.rda')
codega_data <- read.csv(file = 'C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/aNSC_qNSC_profiles.csv', header = T)
#sequencing data
seq_data <- load('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/DATASET_070717.rda')
seq_data <- DATR
row.names(seq_data) <- seq_data$ID
seq_data <- seq_data[,-1]

#sample details
pheno_seq <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/Samples_Davies_040717.txt', stringsAsFactors = F)
row.names(pheno_seq) <- pheno_seq$ID

################test
#this tests if the condition effect (coculture) is different in KO compared to flox##############
pval_threshold <- 0.05
coculture_KO <- results(filtered_davies, name = 'genotypeKO.treatmentendo')
sig_KO <- coculture_KO[coculture_KO$padj < pval_threshold & !is.na(coculture_KO$padj),]
sig_KO <- as.data.frame(sig_KO) #only 36 genes o.o
write.table(sig_KO, 'weird_genes.txt', sep = '\t')



#this tests if the condition effect (BMP) is different in KO compared to flox##############

bmp_KO <- results(filtered_davies, name = 'genotypeKO.treatmentbmp')
sig_bmp <- bmp_KO[bmp_KO$padj < pval_threshold & !is.na(bmp_KO$padj),]
sig_bmp <- as.data.frame(sig_bmp) #only 36 genes o.o
write.table(sig_bmp, 'weird_genes_bmp.txt', sep = '\t')


###############or we can just construct weird factors######################
davies_pdata <- data.frame(row.names = colnames(seq_data), conditions = conditions, genotype = genotype, treatment = treatment, 
													 replicate = as.factor(pheno_seq$Replicate), group = factor(paste0(genotype, treatment)) )

davies_ds_grouped <- DESeqDataSetFromMatrix(countData = seq_data, colData = davies_pdata, design = ~ group)
keep <- rowSums(DESeq2::counts(davies_ds_grouped)) >= 10
filtered_davies_gr <- davies_ds_grouped[keep,] #Remove genes with less that 10 counts in total in all samples
filtered_davies_gr$treatment <- relevel(filtered_davies_gr$group, ref = "floxcontrol_endo")
filtered_davies_gr$genotype <- relevel(filtered_davies_gr$genotype, ref = "flox")

grouped_ds <- DESeq(filtered_davies_gr)
save(grouped_ds,file = 'grouped_DESeq.RData')
pdf('figure2.pdf', height = 25, width = 25)
#####################Coculture WT vs coculture KO#############
KO_cc <- results(grouped_ds, contrast=c('group','KOendo','floxendo'))
sig_KO_cc <- KO_cc[KO_cc$padj < pval_threshold & !is.na(KO_cc$padj) & abs(KO_cc$log2FoldChange) > 0,]
sig_KO_cc <- as.data.frame(sig_KO_cc)

up <- sig_KO_cc[sig_KO_cc$log2FoldChange > 0,]
down <- sig_KO_cc[sig_KO_cc$log2FoldChange < 0 ,]

up_entrez <- bitr(row.names(up), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
down_entrez <- bitr(row.names(down), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)


up_kegg <- enrichKEGG(up_entrez$ENTREZID, org = 'mmu')
down_kegg <- enrichKEGG(down_entrez$ENTREZID, org = 'mmu')
enrichMap(up_kegg)
enrichMap(down_kegg)##fatty acid degradation down regulated? is that what we want?

foxo <- up_kegg[grep(up_kegg$Description, pattern = 'Fox', fixed = T),]
foxo_id <- strsplit(foxo$geneID, '/', fixed = T, perl = F)
foxo_id <- bitr(foxo_id[[1]], fromType = 'ENTREZID',toType = 'GENENAME', OrgDb = org.Mm.eg.db)

write.table(KO_cc, 'KO_coculture_vs_flox_coculture.csv', sep =',')

############################Check overlap with p53 lists##############
p53_targets <- read.csv('../06_Analysis/02_Data/TFCheckpoint_download_180515.csv', header = T, strings = F)
p53_targets$gene_symbol<- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(p53_targets$gene_symbol), perl=TRUE)

p53_up <- intersect(p53_targets$gene_symbol, up_entrez$ALIAS)
p53_down <- intersect(p53_targets$gene_symbol, down_entrez$ALIAS)

p53_up_entrez <- bitr(p53_up, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
p53_down_entrez <- bitr(p53_down, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

p53_up_kegg <- enrichKEGG(p53_up_entrez$ENTREZID, org = 'mmu')
p53_down_kegg <- enrichKEGG(p53_down_entrez$ENTREZID, org = 'mmu')
enrichMap(p53_up_kegg)
enrichMap(p53_down_kegg)
####################oncogene paper#############
onc_p53 <- read.delim('oncogene_p53_targets.txt', header = T, strings = F)
onc_p53$Gene.Symbol<- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(onc_p53$Gene.Symbol), perl=TRUE)

onc_p53_up <- intersect(onc_p53$Gene.symbol, up_entrez$ALIAS)
write.table(onc_p53_up, 'up_regulated_onc_p53.txt', sep = '\t')
onc_p53_down <- intersect(onc_p53$Gene.Symbol, down_entrez$ALIAS)
write.table(onc_p53_down, 'down_regulated_onc_p53.txt', sep = '\t')

onc_p53_down_entrez <- bitr(onc_p53_down, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
kegg <- enrichKEGG(onc_p53_down_entrez$ENTREZID, org = 'mmu')
enrichMap(kegg)
dev.off()
################################Up regulated and downregulated plot#############
up_bub <- i_plot_bubbles(up_kegg, KO_cc, 5, T, 'up') 
down_bub <- i_plot_bubbles(down_kegg, KO_cc, 5, T, 'down') 

all_bub <- join_bubble_objects(up_bub, down_bub)
p <- ggplot(all_bub, aes(x = means, y = -log10(pvalue), size = counter, fill = means, label = term))
p2<- p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel(size = 4, color = 'black') + theme_light() + scale_fill_distiller(palette = 'Spectral')

#############################Clustering data############
vst_seq <- vst(grouped_ds)
vst_seq <- as.data.frame(assay(vst_seq))
only_int <- vst_seq[which(row.names(vst_seq) %in% row.names(sig_KO_cc)),]
log2_int <- log2_fold_change(only_int, samples = c(4,8,12),reference = c(2,6,10)) 
col_pretty <- colorRampPalette(c("blue", "azure1", "red"))
table_cluster <- pheatmap(log2_int, color = col_pretty(15), breaks = seq(-1,1,length.out = 15), clustering_method = 'ward.D2')
cluster_list <- table_cluster$kmeans$cluster

###############################plot genes in stat1#############
stat_targets <- read.csv('../GRSB-7-2013-041-s001/Supplementary Files 11433/stat1_targets.csv', header = T, strings = F, skip =1)
stat_targets$Gene.Symbol<- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(stat_targets$Gene.Symbol), perl=TRUE)

stat_up <- intersect(stat_targets$Gene.Symbol, up_entrez$ALIAS)
stat_down <- intersect(stat_targets$Gene.Symbol, down_entrez$ALIAS)

stat_up_entrez <- bitr(stat_up, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
stat_down_entrez <- bitr(stat_down, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

stat_up_kegg <- enrichKEGG(stat_up_entrez$ENTREZID, org = 'mmu')
stat_down_kegg <- enrichKEGG(stat_down_entrez$ENTREZID, org = 'mmu')
enrichMap(stat_up_kegg)
enrichMap(stat_down_kegg)


#############################Check codega genes in weird genes################
codega_data <- read.csv(file = 'C:/Users/am4613/OneDrive - Imperial College London/Simona_project/7_Analysis/02_Data/aNSC_qNSC_profiles.csv', header = T)
weird_qNSC <- intersect(codega_data$qNSC, row.names(sig_KO))
weird_aNSC <- intersect(codega_data$aNSC, row.names(sig_KO))

test_bub <- i_plot_bubbles(rbind(as.data.frame(up_kegg), as.data.frame(down_kegg)), KO_cc,5,T,'up')

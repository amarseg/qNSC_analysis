##########################Create GO table################
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
source('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/Grouping.R')
source('P:/Simona_project/qNSC_analysis/useful_funs.R')
setwd('P:/Simona_project/')
# load('figure_2/grouped_DESeq.RData')
# #load all data
# coculture_de <- read.csv('figure_1/coculture_deseq.csv', header = T, strings = F)
# p53_coculture <- read.csv('figure_2/KO_coculture_vs_flox_coculture.csv', header = T, strings = F)
# bmp_de <- read.csv('figure_3/bmp_deseq.csv', header = T, strings = F)
# p53KO_de <- read.csv('figure_4/p53KO_de.csv', header = T, strings = F)
# bmpp53KO_de <- read.csv('figure_4/ko_bmp_vs_flox_bmp.csv', header = T, strings = F)

seq_data <- load('P:/Simona_project/06_Analysis/02_Data/DATASET_070717.rda')
seq_data <- DATR[,-1]
row.names(seq_data) <- DATR[,1]
codega_data <- read.csv(file = 'C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/aNSC_qNSC_profiles.csv', header = T)

source('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/Grouping.R')
pheno_seq <- read.delim('P:/Simona_project/06_Analysis/02_Data/Samples_Davies_040717.txt', stringsAsFactors = F)

davies_pdata <- data.frame(row.names = colnames(seq_data), conditions = conditions, genotype = genotype, treatment = treatment, 
                           replicate = as.factor(pheno_seq$Replicate))

davies_ds <- DESeqDataSetFromMatrix(countData = seq_data, colData = davies_pdata, design = ~ replicate)
davies_ds$group <- factor(paste0(davies_ds$genotype, ": ", davies_ds$treatment))
design(davies_ds) <- formula(~genotype + treatment + treatment:genotype) #taken from DESeq2 vignette

keep <- rowSums(DESeq2::counts(davies_ds)) >= 10
filtered_davies <- davies_ds[keep,] #Remove genes with less that 10 counts in total in all samples
filtered_davies$treatment <- relevel(filtered_davies$treatment, ref = "control_endo")
filtered_davies$genotype <- relevel(filtered_davies$genotype, ref = "flox")

filtered_davies <- DESeq(filtered_davies)




davies_pdata <- data.frame(row.names = colnames(seq_data), conditions = conditions, genotype = genotype, treatment = treatment, 
                           replicate = as.factor(pheno_seq$Replicate), group = factor(paste0(genotype, treatment)) )

davies_ds_grouped <- DESeqDataSetFromMatrix(countData = seq_data, colData = davies_pdata, design = ~ group)
keep <- rowSums(DESeq2::counts(davies_ds_grouped)) >= 10
filtered_davies_gr <- davies_ds_grouped[keep,] #Remove genes with less that 10 counts in total in all samples
filtered_davies_gr$treatment <- relevel(filtered_davies_gr$group, ref = "floxcontrol_endo")
filtered_davies_gr$genotype <- relevel(filtered_davies_gr$genotype, ref = "flox")

grouped_ds <- DESeq(filtered_davies_gr)
#extract p53KO coculture vs. p53KO
p53KO_cc_vs_p53KO <- results(grouped_ds, contrast = c('group','KOendo','KOcontrol_endo')) %>%
  as.data.frame()

p53_coculture <- results(grouped_ds, contrast=c('group','KOendo','floxendo')) %>%
  as.data.frame()

bmp_de <- results(grouped_ds, contrast=c('group','floxbmp','floxcontrol_bmp')) %>%
  as.data.frame()

bmpp53KO_de <- results(grouped_ds, contrast=c('group','KObmp','floxbmp')) %>%
  as.data.frame()

p53KO_de <- results(grouped_ds, contrast=c('group','KOcontrol_endo','floxcontrol_endo')) %>%
  as.data.frame()

coculture_de <- results(filtered_davies, contrast = c('treatment','endo','control_endo')) %>%
  as.data.frame()


colnames(coculture_de) <- paste0(colnames(coculture_de), '_coculture_vs_control')
colnames(p53_coculture) <- paste0(colnames(p53_coculture), '_p53KO_coculture_vs_flox_coculture')
colnames(bmp_de) <- paste0(colnames(bmp_de), '_bmp_treated_vs_control')
colnames(p53KO_de) <- paste0(colnames(p53KO_de), '_p53KO_vs_control')
colnames(bmpp53KO_de) <- paste0(colnames(bmpp53KO_de), '_p53KO_bmp_vs_flox_bmp')
colnames(p53KO_cc_vs_p53KO) <- paste0(colnames(p53KO_cc_vs_p53KO), '_p53KO_coculture_vs_p53_control')

coculture_de <- rownames_to_column(coculture_de)
p53_coculture <- rownames_to_column(p53_coculture)
bmp_de <- rownames_to_column(bmp_de)
p53KO_de <- rownames_to_column(p53KO_de)
bmpp53KO_de <- rownames_to_column(bmpp53KO_de)
p53KO_cc_vs_p53KO <- rownames_to_column(p53KO_cc_vs_p53KO)


coculture_de <- coculture_de[,c(1,3,5,7)]
p53_coculture <- p53_coculture[,c(1,3,5,7)]
bmp_de <- bmp_de[,c(1,3,5,7)]
p53KO_de <- p53KO_de[,c(1,3,5,7)]
bmpp53KO_de <- bmpp53KO_de[,c(1,3,5,7)]
p53KO_cc_vs_p53KO <- p53KO_cc_vs_p53KO[,c(1,3,5,7)]


all_DE <- coculture_de %>% 
  inner_join(p53_coculture, by = 'rowname') %>% 
  inner_join(bmp_de,by = 'rowname') %>% 
  inner_join(p53KO_de, by = 'rowname') %>% 
  inner_join(bmpp53KO_de, by = 'rowname') %>% 
  inner_join(p53KO_cc_vs_p53KO, by = 'rowname')

####I need to add the normalised counts now

pheno_seq <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/Samples_Davies_040717.txt', stringsAsFactors = F)
row.names(pheno_seq) <- pheno_seq$ID
sz_filtered <-estimateSizeFactors(filtered_davies_gr)
norm_counts <- as.data.frame(DESeq2::counts(sz_filtered, normalized = T))
norm_counts <- rownames_to_column(norm_counts)
###annotate using codega tables
codega_data <- read.csv(file = 'C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/aNSC_qNSC_profiles.csv', header = T)

codega_data[codega_data == ""] <- NA
m_codega <- gather(codega_data,na.rm = T)
#I remove duplicates, how can this people upload a gene set with duplicates?!
names_split <- strsplit(m_codega$value, split = ' /// ', fixed = T)
good_names <- unlist(lapply(names_split, function(x){x[1]}))
m_codega$value <- good_names
m_codega <- m_codega[!duplicated(m_codega$value),]

row.names(m_codega) <- m_codega$value
m_codega <- dplyr::select(m_codega, key)
m_codega <- rownames_to_column(m_codega)

all_de_codega <- left_join(all_DE, m_codega, by = 'rowname')
###add bobadilla o boadilla
boa_data <- read.delim(file = 'C:/Users/am4613/OneDrive - Imperial College London/Simona_project/Boadilla_supp/bobadilla_clusters.txt', header = T, stringsAsFactors = F)
#sequencing data
aNSC <- boa_data[which(boa_data$cluster == 4 | boa_data$cluster == 5 | boa_data$cluster == 6),]$gene_symbol
qNSC <- boa_data[which(boa_data$cluster == 2 | boa_data$cluster == 3),]$gene_symbol

boadilla <- data.frame(genes = c(aNSC, qNSC), key = c(rep('aNSC',1026), rep('qNSC',315)))
boadilla$genes<- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(boadilla$genes), perl=TRUE)
m_boadilla <- boadilla[which(boadilla$genes != ''),]
row.names(m_boadilla) <- m_boadilla$genes

all_de_boa <- left_join(all_de_codega, m_boadilla, by = c('rowname' = 'genes'))

####last but not least, add the normalised counts
all_de_codega_counts <- left_join(all_de_boa, norm_counts, by = 'rowname')
colnames(all_de_codega_counts)[1] <- 'Gene_names'
colnames(all_de_codega_counts)[c(20,21)] <- c('Codega et al.','Boadilla et al.')
#write.table(all_de_codega_counts, 'supp1.txt', sep = '\t', row.names = F)
##############################Same with enrichments#########################
enricher_function <- function(de_obj, comparison)
{
	require('clusterProfiler')
	pval_thr <- 0.05
	sig <- de_obj[de_obj$padj < pval_thr & !is.na(de_obj$padj),]
	up <- sig[sig$log2FoldChange > 0 ,]
	down <- sig[sig$log2FoldChange < 0 ,]
	
	up_entrez <- bitr(up$rowname, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
	down_entrez <- bitr(down$rowname, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
	
	up_GO <- enrichGO(up_entrez$ENTREZID, OrgDb = org.Mm.eg.db, keyType = 'ENTREZID')
	up_GO <- as.data.frame(up_GO)
	up_GO$database <- 'GO'
	up_GO$direction <- 'Upregulated'
	down_GO <- enrichGO(down_entrez$ENTREZID, OrgDb = org.Mm.eg.db, keyType = 'ENTREZID')
	down_GO <- as.data.frame(down_GO)
	down_GO$database <- 'GO'
	down_GO$direction <- 'Downregulated'
	
	up_david <- enrichDAVID(up_entrez$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk')
	up_david <- as.data.frame(up_david)
	up_david$database <- 'DAVID'
	up_david$direction <- 'Upregulated'
	down_david <- enrichDAVID(down_entrez$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk')
	down_david <- as.data.frame(down_david)
	down_david$database <- 'DAVID'
	down_david$direction <- 'Downregulated'
	
	up_KEGG <- enrichKEGG(up_entrez$ENTREZID, organism = 'mmu')
	up_KEGG <- as.data.frame(up_KEGG)
	up_KEGG$database <- 'KEGG'
	up_KEGG$direction <- 'Upregulated'
	down_KEGG <- enrichKEGG(down_entrez$ENTREZID, organism = 'mmu')
	down_KEGG <- as.data.frame(down_KEGG)
	down_KEGG$database <- 'KEGG'
	down_KEGG$direction <- 'Downregulated'
	
	all_enrich <- rbind(up_GO, down_GO, down_KEGG, up_KEGG)	
	all_enrich$Comparison <- comparison
	
	split_ids <- strsplit(all_enrich$geneID, split = '/', fixed = T)
	split_symbols <- lapply(split_ids, FUN = bitr, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = org.Mm.eg.db)
	split_symb <- lapply(split_symbols, function(x){x$SYMBOL})
	col_symb <- as.vector(lapply(split_symb, FUN = paste, collapse = '/'))
	all_enrich$geneID <- unlist(col_symb)
	return(all_enrich)
}

coculture_enrich <- enricher_function(coculture_de, comparison = 'coculture_vs_control')
p53coculture_enrich <- enricher_function(p53_coculture, comparison = 'p53_KOcoculture_vs_flox_coculture')
bmp_enrich <- enricher_function(bmp_de, comparison = 'bmp_vs_control')
p53KO_enrich <- enricher_function(p53KO_de, comparison = 'p53KO_vs_flox')
bmpp53KO_enrich <- enricher_function(bmpp53KO_de, comparison = 'p53_KO_bmp_vs_flox_bmp')
p53_KO_cc_p53KO_enrich <- enricher_function(p53KO_cc_vs_p53KO, comparison = '_p53KO_coculture_vs_p53_control')

all_enrichments <- rbind(coculture_enrich, p53coculture_enrich, bmp_enrich, p53KO_enrich, bmpp53KO_enrich, p53_KO_cc_p53KO_enrich )

#write.table(all_enrichments, 'enrichments_table.txt', sep = '\t', row.names = F)

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
setwd('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/figure_1_bobadilla')
###############Load data######
#Load data
load('../DESeq_filtered.rda')
codega_data <- read.delim(file = 'C:/Users/am4613/OneDrive - Imperial College London/Simona_project/Boadilla_supp/bobadilla_clusters.txt', header = T, stringsAsFactors = F)
#sequencing data
seq_data <- load('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/DATASET_070717.rda')
seq_data <- DATR
row.names(seq_data) <- seq_data$ID
seq_data <- seq_data[,-1]
#sample details
pheno_seq <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/Samples_Davies_040717.txt', stringsAsFactors = F, skip = 1)
row.names(pheno_seq) <- pheno_seq$ID

###Get DE genes from results object
p_adj_thr <- 0.05
res_cc <- results(filtered_davies, filterFun = ihw, contrast = c('treatment','endo','control_endo'))
cc_genes <- subset(res_cc, padj < p_adj_thr)
cc_sig_lfc <- subset(res_cc, padj < p_adj_thr & abs(log2FoldChange) > 0)

###Calculate fold changes for all three replicates
interesting <- c(1,2,5,6,9,10)
davies_rlog <- vst(filtered_davies)
davies_vst <- as.data.frame(assay(davies_rlog))
only_cc_seq <- davies_vst[, interesting]

fsc_cc <- davies_vst[,c(2,4,6)]/davies_vst[,c(1,3,5)]
fsc_cc <- fsc_cc +0.00001
log2_fsc <- log2(fsc_cc)
#fold change heatmap
col_pretty <- colorRampPalette(c("blue", "azure1", "red"))
pheatmap(log2_fsc, show_rownames = F, color = col_pretty(30))##all_genes!!!

##create table for annotation
aNSC <- codega_data[which(codega_data$cluster == 4 | codega_data$cluster == 5 | codega_data$cluster == 6),]$gene_symbol
qNSC <- codega_data[which(codega_data$cluster == 2 | codega_data$cluster == 3),]$gene_symbol

boadilla <- data.frame(genes = c(aNSC, qNSC), key = c(rep('aNSC',1026), rep('qNSC',315)))
boadilla$genes<- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(boadilla$genes), perl=TRUE)
m_boadilla <- boadilla[which(boadilla$genes != ''),]
row.names(m_boadilla) <- m_boadilla$genes
m_boadilla <- dplyr::select(m_boadilla, key)

##################Heatmap#####
#Heatmap code
##only DE genes
fsc_cc_DE <- log2_fsc[which(row.names(log2_fsc) %in% row.names(cc_sig_lfc)),]
pheatmap(fsc_cc_DE, show_rownames = F, color = col_pretty(30), annotation_row = m_boadilla, show_colnames = F, clustering_method = 'ward.D2', filename = 'iamaplot.pdf')##onlyDE genes

test <- fsc_cc_DE[which(row.names(fsc_cc_DE) %in% row.names(m_codega)),]
pheatmap(test, show_rownames = F, color = col_pretty(30), annotation_row = m_boadilla, show_colnames = F, clustering_method = 'ward.D2', filename = 'iamaplot.pdf')##onlyDE genes

########################Barplot#####
#Bar plots showing what goes in which direction

up_reg <- subset(cc_sig_lfc, log2FoldChange >0 )
down_reg <- subset(cc_sig_lfc, log2FoldChange < 0 )

m_codega$direction <- 'NA'
m_codega[which(row.names(m_codega) %in% row.names(up_reg) & m_codega$key == 'qNSC'),]$direction <- 'Right Direction'
m_codega[which(row.names(m_codega) %in% row.names(down_reg) & m_codega$key == 'aNSC'),]$direction <- 'Right Direction'

m_codega[which(row.names(m_codega) %in% row.names(down_reg) & m_codega$key == 'qNSC'),]$direction <- 'Wrong Direction'
m_codega[which(row.names(m_codega) %in% row.names(up_reg) & m_codega$key == 'aNSC'),]$direction <- 'Wrong Direction'

m_codega[which(m_codega$direction == 'NA'),]$direction <- 'Not Differentially regulated'

m_codega$direction <- factor(m_codega$direction, levels = c('Wrong Direction','Not Differentially regulated','Right Direction'))

write.table(m_codega, 'regulation_codega.txt', sep = '\t')

p <- ggplot(m_codega, aes(key, fill = direction))
p2 <- p + geom_bar(position = "fill") + scale_fill_brewer(type = 'qual', palette = 'Dark2') + theme_light() 
ggsave('barplot_DE.pdf',plot = p2, device = 'pdf')



####################up and down barplot######
codega_up_down <- m_codega
codega_up_down$direction <- 'NA'
codega_up_down[which(row.names(codega_up_down) %in% row.names(up_reg) & codega_up_down$key == 'qNSC'),]$direction <- 'Up'
codega_up_down[which(row.names(codega_up_down) %in% row.names(down_reg) & codega_up_down$key == 'aNSC'),]$direction <- 'Down'

codega_up_down[which(row.names(codega_up_down) %in% row.names(down_reg) & codega_up_down$key == 'qNSC'),]$direction <- 'Down'
codega_up_down[which(row.names(codega_up_down) %in% row.names(up_reg) & codega_up_down$key == 'aNSC'),]$direction <- 'Up'

codega_up_down[which(codega_up_down$direction == 'NA'),]$direction <- 'Not Differentially regulated'

codega_up_down$direction <- factor(codega_up_down$direction, levels = c('Up','Down','Not Differentially regulated'))

p <- ggplot(codega_up_down, aes(key, fill = direction))
p2 <- p + geom_bar(position = "fill") + scale_fill_brewer(type = 'qual', palette = 'Dark2') + theme_light() 
ggsave('barplot_DE_up_down.pdf',plot = p2, device = 'pdf')

#########################Waste of time v1.0#####
#I am going to try to make a piechart or another barplot of GO terms
m_codega$up_down <- codega_up_down$direction

comparing_codega <- rownames_to_column(m_codega)

entrez_id <- bitr(comparing_codega$rowname, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

codega_entrez <- inner_join(comparing_codega, entrez_id, by = c('rowname' = 'ALIAS'))
# 
# formula_res <- _w
# dp<-dotplot(formula_res, x=~direction, showCategory = 15 ) + ggplot2::facet_grid(~key)
# ggsave('dotplot_aNSC_qNSC.pdf', plot = dp, device = 'pdf', scale = 3.25)
# 
# ##Calulate similarities between GO categories
# library('GOSemSim')
# library('ggrepel')
# similarity_mmu <- godata(OrgDb = org.Mm.eg.db, keytype = "ENTREZID", ont='MF', computeIC = TRUE)
# simp_res <- clusterProfiler::simplify(formula_res, cutoff = 0.7, by = "p.adjust",select_fun = min, measure = "Rel", semData = similarity_mmu)
# dp <- dotplot(simp_res, x=~direction, showCategory = 15) + ggplot2::facet_grid(~key)
# ggsave('dotplot_simplified.pdf', plot = dp, device = 'pdf', scale = 3.25)
# 
# ################
# #This is how you wasted 3 hours trying to make a plot
# df <- as.data.frame(formula_res)
# df$plot.value <- -log10(df$p.adjust)
# p <- ggplot(df, aes(x = direction, y = plot.value, size = Count, label = Description))
# p2<- p + geom_point(shape = 21, colour = "#000000", fill = "#40b8d0", alpha = 0.75) 
# p2 + facet_grid(~key) +  geom_text_repel(data=subset(df, plot.value > 3), size = 3, alpha = 0.5)
# 
# p <- ggplot(df, aes(x = direction, y = Count, size = plot.value, label = Description))
# p2<- p + geom_point(shape = 21, colour = "#000000", fill = "#40b8d0", alpha = 0.75) 
# p2 + facet_grid(~key) +  geom_text_repel(data=subset(df, plot.value > 3), size = 4, alpha = 0.5)

############Waste of time v2.0####
#Pretty GO plots Colours need to be changed!!!!!
library('ggrepel')
right <- codega_entrez[codega_entrez$key == 'qNSC' & codega_entrez$direction == 'Right Direction',]
right_dv <- enrichDAVID(right$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk')
enrichMap(right_dv)

wrong <- codega_entrez[codega_entrez$key == 'qNSC' & codega_entrez$direction == 'Not Differentially regulated',]
wrong_dv <- enrichDAVID(wrong$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk')
enrichMap(wrong_dv)

dv_df_r <- as.data.frame(right_dv)
dv_df_r$plot.value <- -log10(dv_df_r$p.adjust)
p <- ggplot(dv_df_r, aes(x = qvalue, y = plot.value, size = Count, label = Description))
p2<- p + geom_point(shape = 21, colour = "#000000", fill = "#40b8d0", alpha = 0.75) 
p2 +  geom_text_repel( data=subset(dv_df_r, p.adjust < 0.0002),size = 4, alpha = 0.5)

dv_df_w <- as.data.frame(wrong_dv)
dv_df_w$plot.value <- -log10(dv_df_w$p.adjust)
p <- ggplot(dv_df_w, aes(x = qvalue, y = plot.value, size = Count, label = Description))
p2<- p + geom_point(shape = 21, colour = "#000000", fill = "#40b8d0", alpha = 0.75) 
p2 +  geom_text_repel( data=subset(dv_df_w, p.adjust < 0.0002),size = 4, alpha = 0.5)

r_w_list <- list(right$ENTREZID, wrong$ENTREZID)
names(r_w_list) <- c('right','Wrong')
m_kegg_comp <- compareCluster(r_w_list, 'enrichKEGG', organism = 'mmu')
dotplot(m_kegg_comp)

kegg_df <- as.data.frame(m_kegg_comp)
kegg_df$plot.value <- -log10(kegg_df$p.adjust)
p <- ggplot(kegg_df, aes(x = plot.value, y = Cluster, size = Count, label = Description, fill = Cluster))
p2<- p + geom_point(shape = 21,  alpha = 0.75) 
p2 +  geom_text_repel(data= subset(kegg_df, p.adjust < 0.01),size = 4) + scale_fill_brewer(palette = 'Dark2') + theme_light()
ggsave('qNSC_enrichment.pdf')


right_a <- codega_entrez[codega_entrez$key == 'aNSC' & codega_entrez$direction == 'Right Direction',]
wrong_a <- codega_entrez[codega_entrez$key == 'aNSC' & codega_entrez$direction == 'Not Differentially regulated',]
r_w_list_a <- list(right_a$ENTREZID, wrong_a$ENTREZID)
names(r_w_list_a) <- c('right','Wrong')
m_kegg_comp_a <- compareCluster(r_w_list_a, 'enrichKEGG', organism = 'mmu')
dotplot(m_kegg_comp_a)

kegg_df_a <- as.data.frame(m_kegg_comp_a)
kegg_df_a$plot.value <- -log10(kegg_df_a$p.adjust)
p <- ggplot(kegg_df_a, aes(x = plot.value, y = Cluster, size = Count, label = Description, fill = Cluster))
p2<- p + geom_point(shape = 21,  alpha = 0.75) 
p2 +  geom_text_repel(data= subset(kegg_df_a, p.adjust < 0.01),size = 4) + scale_fill_brewer(palette = 'Dark2') + theme_light()
ggsave('aNSC_enrichment.pdf')


#######GOplot#####
library(GOplot)
#Need to create things for circle_dat
#terms = r 'category', 'ID', 'term', adjusted p-value ('adj_pval') and 'genes'
term_test <- as.data.frame(right_dv)
term_test$geneID <- gsub(term_test$geneID, pattern = '/', replacement = ', ', fixed = T )

term_thing <- data.frame(category = 'BP', ID = term_test$ID, adj_pval = term_test$p.adjust, genes = term_test$geneID,
												 term = term_test$Description)


#genes 'ID', 'logFC'

genes_go <- data.frame(logFC=cc_sig_lfc$log2FoldChange, SYMBOL = row.names(cc_sig_lfc)) 
entrez_ids <- bitr(row.names(cc_sig_lfc), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db )

all_genes <- inner_join(genes_go, entrez_ids, by = c('SYMBOL' = 'ALIAS'))
colnames(all_genes)[3] <- 'ID'
###combine the two?

circ <- circle_dat(terms = term_thing, genes = all_genes)
GOBubble(circ, labels = 3, ID = F)

############GOplot function############
#The function takes an enrichment object from clusterProfiler
#and a results object from DESeq2
i_plot_bubbles <- function(enrichment_obj, results_obj, labels, ID, title)
{
	#Prepare enrichment object for circle_dat function
	term_test <- as.data.frame(enrichment_obj)
	term_test$geneID <- gsub(term_test$geneID, pattern = '/', replacement = ', ', fixed = T )
	
	term_thing <- data.frame(category = 'BP', ID = term_test$ID, adj_pval = term_test$p.adjust, genes = term_test$geneID,
													 term = term_test$Description)
	
	genes_go <- data.frame(logFC=results_obj$log2FoldChange, SYMBOL = row.names(results_obj)) 
	entrez_ids <- bitr(row.names(results_obj), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db )
	
	all_genes <- inner_join(genes_go, entrez_ids, by = c('SYMBOL' = 'ALIAS'))
	colnames(all_genes)[3] <- 'ID'
	###combine the two?
	
	circ <- circle_dat(terms = term_thing, genes = all_genes)
	GOBubble(circ, labels = labels, ID = ID, title = title)
	return(circ)
}
kegg_df$Cluster <- as.factor(kegg_df$Cluster)
qNSC_kegg <- base::split(kegg_df, f = kegg_df$Cluster)
pdf('bubble_plots_kegg.pdf', width = 14, height = 14)
i_plot_bubbles(qNSC_kegg[[1]], res_cc, 3, F, title = 'qNSC_right')
t <- i_plot_bubbles(qNSC_kegg[[2]], res_cc, 3, F, title = 'qNSC_wrong')

aNSC_kegg <- base::split(kegg_df_a, f = kegg_df_a$Cluster)
i_plot_bubbles(aNSC_kegg[[1]],res_cc , 3, F, title = 'aNSC_righ')
i_plot_bubbles(aNSC_kegg[[2]],res_cc , 3, F, title = 'aNSC_righ')

dev.off()
#i_plot_bubbles(aNSC_kegg[[2]], res_cc, 3, F)

m_david_comp <- compareCluster(r_w_list, 'enrichDAVID', david.user = 'a.martinez-segura13@imperial.ac.uk')
m_david_comp_a <- compareCluster(r_w_list_a, 'enrichDAVID', david.user = 'a.martinez-segura13@imperial.ac.uk')

m_david_comp <- as.data.frame(m_david_comp)
m_david_comp_a <- as.data.frame(m_david_comp_a)

qNSC_david <- base::split(m_david_comp, f = m_david_comp$Cluster)
pdf('bubble_plots_DAVID.pdf', width = 14, height = 14)
#i_plot_bubbles(qNSC_david[[1]], res_cc, 3, F, title = 'qNSC_right')
t <- i_plot_bubbles(qNSC_david[[2]], res_cc, 3, F, title = 'qNSC_wrong')

aNSC_david <- base::split(m_david_comp_a, f = m_david_comp_a$Cluster)
i_plot_bubbles(aNSC_david[[1]],res_cc , 15, F, title = 'aNSC_right')
i_plot_bubbles(aNSC_david[[2]], res_cc, 3, F, title = 'aNSC_wrong')
dev.off()



m_GO_comp <- compareCluster(r_w_list, 'enrichGO', OrgDb = org.Mm.eg.db, ont = 'BP')
m_GO_comp_a <- compareCluster(r_w_list_a, 'enrichGO', OrgDb = org.Mm.eg.db, ont = 'BP')

m_GO_comp <- as.data.frame(m_GO_comp)
m_GO_comp_a <- as.data.frame(m_GO_comp_a)

qNSC_GO <- base::split(m_GO_comp, f = m_GO_comp$Cluster)
pdf('bubble_plots_GO.pdf', width = 14, height = 14)
i_plot_bubbles(qNSC_GO[[1]], res_cc, 5, F, title = 'qNSC_right')
t <- i_plot_bubbles(qNSC_GO[[2]], res_cc, 5, F, title = 'qNSC_wrong')

aNSC_GO <- base::split(m_GO_comp_a, f = m_GO_comp_a$Cluster)
#i_plot_bubbles(aNSC_GO[[1]],res_cc , 5, F, title = 'aNSC_right')
#i_plot_bubbles(aNSC_GO[[2]], res_cc, 5, F, title = 'aNSC_wrong')
dev.off()




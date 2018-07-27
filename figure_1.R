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
library('pathview')
setwd('P:/Simona_project/figure_1')
###############Load data######
#Load data
load('../DESeq_filtered.rda')
codega_data <- read.csv(file = 'P:/Simona_project/06_Analysis/02_Data/aNSC_qNSC_profiles.csv', header = T)
#sequencing data
seq_data <- load('P:/Simona_project/06_Analysis/02_Data/DATASET_070717.rda')
seq_data <- DATR
row.names(seq_data) <- seq_data$ID
seq_data <- seq_data[,-1]
#sample details
pheno_seq <- read.delim('P:/Simona_project/06_Analysis/02_Data/Samples_Davies_040717.txt', stringsAsFactors = F)
row.names(pheno_seq) <- pheno_seq$ID

###Get DE genes from results object
p_adj_thr <- 0.05
res_cc <- results(filtered_davies, filterFun = ihw, contrast = c('treatment','endo','control_endo'))
cc_genes <- subset(res_cc, padj < p_adj_thr)
cc_sig_lfc <- subset(res_cc, padj < p_adj_thr & abs(log2FoldChange) > 0)
write.table(res_cc, 'coculture_deseq.csv', sep = ',')
###Calculate fold changes for all three replicates
interesting <- c(1,2,5,6,9,10)
davies_vst <- as.data.frame(DESeq2::counts(filtered_davies, normalized = T))
davies_vst <- davies_vst + 0.0000001
fsc_cc <- davies_vst[,c(2,6,10)]/davies_vst[,c(1,5,9)]
log2_fsc <- log2(fsc_cc)
#fold change heatmap
col_pretty <- colorRampPalette(c("blue", "azure1", "red"))
pheatmap(log2_fsc, show_rownames = F, color = col_pretty(30))##all_genes!!!

##create table for annotation
codega_data[codega_data == ""] <- NA
m_codega <- gather(codega_data,na.rm = T)
#I remove duplicates, how can this people upload a gene set with duplicates?!
names_split <- strsplit(m_codega$value, split = ' /// ', fixed = T)
good_names <- unlist(lapply(names_split, function(x){x[1]}))
m_codega$value <- good_names
m_codega <- m_codega[!duplicated(m_codega$value),]

row.names(m_codega) <- m_codega$value
m_codega <- dplyr::select(m_codega, key)

##################Heatmap#####
#Heatmap code
##only DE genes
annot_col <- list(CellType = c(qNSC = '#1B9E77', aNSC = "#D95F02"))
fsc_cc_DE <- log2_fsc[which(row.names(log2_fsc) %in% row.names(cc_sig_lfc)),]
pdf('heatmap_DE_codega.pdf')
pheatmap(fsc_cc_DE, show_rownames = F, color = col_pretty(30), annotation_row = m_codega, 
				 annotation_colors = annot_col[1],show_colnames = F, clustering_method = 'ward.D2', breaks = seq(-15,15,by = 1))##onlyDE genes
dev.off()
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
write.table(codega_up_down, 'regulation_codega.txt', sep = '\t')

p <- ggplot(codega_up_down, aes(key, fill = direction))
p2 <- p + geom_bar(position = "fill") + scale_fill_brewer(type = 'qual', palette = 'Dark2') + theme_light() 
ggsave('barplot_DE_up_down.pdf',plot = p2, device = 'pdf')

m_codega$up_down <- codega_up_down$direction
#########################Waste of time v1.0#####
#I am going to try to make a piechart or another barplot of GO terms


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
right_dv <- enrichDAVID(right$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk',annotation = 'GOTERM_BP_3')
enrichMap(right_dv)

wrong <- codega_entrez[codega_entrez$key == 'qNSC' & codega_entrez$direction == 'Not Differentially regulated',]
wrong_dv <- enrichDAVID(wrong$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk', annotation = 'GOTERM_BP_3')
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
write.table(m_kegg_comp, 'qNSC_kegg.txt', sep = '\t')

m_david_comp <- compareCluster(r_w_list, 'enrichDAVID', david.user = 'a.martinez-segura13@imperial.ac.uk', annotation = 'GOTERM_BP_3')
write.table(m_david_comp, 'qNSC_david.txt',sep = '\t')

m_GO_comp <- compareCluster(r_w_list, 'enrichGO', OrgDb = org.Mm.eg.db , ont = 'BP')
write.table(m_GO_comp, 'qNSC_GO.txt',sep = '\t')





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
write.table(m_kegg_comp_a, 'aNSC_kegg.txt',sep = '\t')

m_david_comp_a <- compareCluster(r_w_list_a, 'enrichDAVID', david.user = 'a.martinez-segura13@imperial.ac.uk', annotation = 'GOTERM_BP_3')
write.table(m_david_comp_a, 'aNSC_david.txt',sep = '\t')

m_GO_comp_a <- compareCluster(r_w_list_a, 'enrichGO', OrgDb = org.Mm.eg.db , ont = 'BP')
write.table(m_GO_comp_a, 'aNSC_GO.txt',sep = '\t')


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
right_q_kegg <-i_plot_bubbles(qNSC_kegg[[1]], res_cc, 3, F, title = 'qNSC_right')
write.table(right_q_kegg, 'Up_qNSC.txt', sep = '\t')
wrong_q_kegg <- i_plot_bubbles(qNSC_kegg[[2]], res_cc, 1.5, F, title = 'qNSC_wrong')
write.table(wrong_q_kegg, 'Not_DE_qNSC.txt', sep = '\t')

aNSC_kegg <- base::split(kegg_df_a, f = kegg_df_a$Cluster)
right_a_kegg <- i_plot_bubbles(aNSC_kegg[[1]],res_cc , 3, F, title = 'aNSC_righ')
write.table(right_q_kegg, 'Down_aNSC.txt', sep = '\t')
wrong_a_kegg <- i_plot_bubbles(aNSC_kegg[[2]], res_cc, 3, F)
write.table(aNSC_kegg[[2]], 'Not_DE_aNSC.txt', sep = '\t')
dev.off()


m_david_comp <- compareCluster(r_w_list, 'enrichDAVID', david.user = 'a.martinez-segura13@imperial.ac.uk', annotation = 'GOTERM_BP_3')
m_david_comp_a <- compareCluster(r_w_list_a, 'enrichDAVID', david.user = 'a.martinez-segura13@imperial.ac.uk',annotation = 'GOTERM_BP_3')

m_david_comp <- as.data.frame(m_david_comp)
m_david_comp_a <- as.data.frame(m_david_comp_a)

qNSC_david <- base::split(m_david_comp, f = m_david_comp$Cluster)
pdf('bubble_plots_DAVID.pdf', width = 14, height = 14)
q_david_r <- i_plot_bubbles(qNSC_david[[1]], res_cc, 1, F, title = 'qNSC_right')
q_david_w <- i_plot_bubbles(qNSC_david[[2]], res_cc, 1, F, title = 'qNSC_wrong')

aNSC_david <- base::split(m_david_comp_a, f = m_david_comp_a$Cluster)
a_david_r <- i_plot_bubbles(aNSC_david[[1]],res_cc , 1, F, title = 'aNSC_right')
a_david_w <- i_plot_bubbles(aNSC_david[[2]], res_cc, 1, F, title = 'aNSC_wrong')
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
i_plot_bubbles(aNSC_GO[[1]],res_cc , 5, F, title = 'aNSC_right')
#i_plot_bubbles(aNSC_GO[[2]], res_cc, 5, F, title = 'aNSC_wrong')
dev.off()


#####################Stem cell quiescence signature#############
sig <- read.delim('../rando_table.txt', header = T, strings = F)
sig$Gene <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(sig$Gene), perl=TRUE)

row.names(sig) <- sig$Gene
sig <- dplyr::select(sig, Regulation, Function)

only_sig <- fsc_cc[which(row.names(fsc_cc) %in% row.names(sig)),]

pdf('quiescent_stem_cell_signature.pdf', paper = 'a4r')
pheatmap(only_sig, color = col_pretty(30), annotation_row = sig, clustering_method = 'ward.D2')
dev.off()

##########################PPAR-gamma plot
ppar_gamma <- read.delim('../PPAR_datasets/ppar_gamma_targets_gendev.txt', header =F,strings = F, row.names = 1)

ppar_genes <- log2_fsc[which(row.names(log2_fsc) %in% row.names(ppar_gamma)),]
ppar_gamma$DE <- 'Not DE'
ppar_gamma[which(row.names(ppar_gamma) %in% row.names(cc_genes)),]$DE <- 'DE'

pheatmap(ppar_genes, color = col_pretty(30), annotation_row = ppar_gamma, clustering_method = 'ward.D2')

############################Similarly, i take the ppar alpha ones
ppar_alpha <- read.delim('../PPAR_datasets/ppar_alpha_table.txt', header =F,strings = F, row.names = 1)
ppar_alpha_genes <- log2_fsc[which(row.names(log2_fsc) %in% row.names(ppar_alpha)),]
ppar_alpha$DE <- 'Not DE'
ppar_alpha[which(row.names(ppar_alpha) %in% row.names(cc_genes)),]$DE <- 'DE'

pheatmap(ppar_alpha_genes, color = col_pretty(30), annotation_row = ppar_alpha, clustering_method = 'ward.D2')

####################################ppar delta_beta bitches

ppar_beta <- read.delim('../PPAR_datasets/ppar_beta_targets.txt', header = T, strings =F)
ppar_beta$Gene.Name..manually.assigned.<- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(ppar_beta$Gene.Name..manually.assigned.), perl=TRUE)

ppar_beta_df <- as.data.frame(unique(ppar_beta$Gene.Name..manually.assigned.))
colnames(ppar_beta_df) <- 'DE'
row.names(ppar_beta_df) <- ppar_beta_df[,1]
ppar_beta_df[,1] <- 'Not DE'
ppar_beta_df[which(row.names(ppar_beta_df) %in% row.names(cc_genes)),1] <- 'DE'

ppar_beta_genes <- log2_fsc[which(row.names(log2_fsc) %in% row.names(ppar_beta_df)),]


pheatmap(ppar_beta_genes, color = col_pretty(30), annotation_row = ppar_beta_df, clustering_method = 'ward.D2')

########################Do bubble plot with up and down regulated stuff#####
###I take the result of circle_dat and do my own thing because that's how I roll
up_reg_entrez <- bitr(row.names(up_reg), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
up_reg_kegg <- enrichKEGG(up_reg_entrez$ENTREZID, organism = 'mmu')

down_reg_entrez <- bitr(row.names(down_reg), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
down_reg_kegg <- enrichKEGG(down_reg_entrez$ENTREZID, organism = 'mmu')

down_bub <- i_plot_bubbles(down_reg_kegg, res_cc,5, F, title = 'down_regulated_KEGG') 
write.table(down_bub, 'Downregulated_kegg.txt', sep = '\t')
up_bub <- i_plot_bubbles(up_reg_kegg, res_cc,5, F, title = 'up_regulated_KEGG') 
write.table(up_bub, 'Upregulated_kegg.txt', sep = '\t')

all_bub <- bind_rows(down_bub, up_bub) #Average log change for plotting

all_bub %>% group_by(term) %>% summarise(means = mean(logFC), pvalue = mean(adj_pval), counter = n()) %>% as.data.frame()-> medias

p <- ggplot(medias, aes(x = means, y = -log10(pvalue), size = counter, fill = means, label = term))
p2<- p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel(data=subset(medias, -log10(pvalue) >6), size = 4, color = 'black') + theme_light() + scale_fill_distiller(palette = 'Spectral') + xlim(-2,2)
ggsave('UpDown_KEGG.pdf',p2,scale = 1.5, device = 'pdf')


up_reg_david <- enrichDAVID(up_reg_entrez$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk',annotation = 'GOTERM_BP_3')
down_reg_david <- enrichDAVID(down_reg_entrez$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk', annotation = 'GOTERM_BP_3')

down_bub <- i_plot_bubbles(down_reg_david, res_cc,5, F, title = 'down_regulated_DAVID') 
write.table(down_bub, 'Downregulated_david.txt', sep = '\t')
up_bub <- i_plot_bubbles(up_reg_david, res_cc,5, F, title = 'up_regulated_DAVID') 
write.table(up_bub, 'upregulated_david.txt', sep = '\t')

all_bub <- bind_rows(down_bub, up_bub) #Average log change for plotting

all_bub %>% group_by(term) %>% summarise(means = mean(logFC), pvalue = mean(adj_pval), counter = n()) %>% as.data.frame()-> medias

p <- ggplot(medias, aes(x = means, y = -log10(pvalue), size = counter, fill = means, label = term))
p2<- p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel(data=subset(medias, -log10(pvalue) >5), size = 4, color = 'black') + theme_light() + scale_fill_distiller(palette = 'Spectral')
ggsave('UpDown_david.pdf',p2, device = 'pdf')

#############################Heatmap codega and boadilla#################
boadilla_data <- read.delim(file = 'C:/Users/am4613/OneDrive - Imperial College London/Simona_project/Boadilla_supp/bobadilla_clusters.txt', header = T, stringsAsFactors = F)
aNSC <- boadilla_data[which(boadilla_data$cluster == 4 | boadilla_data$cluster == 5 | boadilla_data$cluster == 6),]$gene_symbol
qNSC <- boadilla_data[which(boadilla_data$cluster == 2 | boadilla_data$cluster == 3),]$gene_symbol

boadilla <- data.frame(genes = c(aNSC, qNSC), key = c(rep('aNSC',1026), rep('qNSC',315)))
boadilla$genes<- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(boadilla$genes), perl=TRUE)
m_boadilla <- boadilla[which(boadilla$genes != ''),]
row.names(m_boadilla) <- m_boadilla$genes
m_boadilla <- dplyr::select(m_boadilla, key)

all_id <- merge(m_codega, m_boadilla, by = 'row.names', all = T)
row.names(all_id) <- all_id$Row.names
m_all_id <- dplyr::select(all_id, key.x, key.y)
colnames(m_all_id) <- c('Codega','Boadilla')

pheatmap(fsc_cc_DE, show_rownames = F, color = col_pretty(30), annotation_row = m_all_id, 
				 annotation_colors = annot_col[1],show_colnames = F, clustering_method = 'ward.D2')

#############################DO aNSC_kegg[[nicely]]
term_test <- as.data.frame(aNSC_kegg[[2]])
term_test$geneID <- gsub(term_test$geneID, pattern = '/', replacement = ', ', fixed = T )

term_thing <- data.frame(category = 'BP', ID = term_test$ID, adj_pval = term_test$p.adjust, genes = term_test$geneID,
												 term = term_test$Description)


#genes 'ID', 'logFC'

genes_go <- data.frame(logFC=res_cc$log2FoldChange, SYMBOL = row.names(res_cc)) 
entrez_ids <- bitr(row.names(res_cc), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db )

all_genes <- inner_join(genes_go, entrez_ids, by = c('SYMBOL' = 'ALIAS'))
colnames(all_genes)[3] <- 'ID'
###combine the two?

aNSC_w <- circle_dat(terms = term_thing, genes = all_genes)



##############################Test all_buble
dark_col <- brewer.pal(3, name = 'Dark2')

all_bub_rw <- bind_rows(aNSC_w, right_a_kegg, right_q_kegg, wrong_q_kegg)
all_bub_rw$type <- c(rep('aNSC Not DE',nrow(aNSC_w)), rep('aNSC Downregulated', nrow(right_a_kegg)),
									rep('qNSC Upregulated', nrow(right_q_kegg)), rep('qNSC Not DE', nrow(wrong_q_kegg)))

all_bub_rw %>% group_by(term,type) %>% summarise(means = mean(logFC, na.rm =T), pvalue = mean(adj_pval, na.rm = T), counter = n()) %>% as.data.frame()-> medias
p <- ggplot(medias, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel(data=subset(medias, -log10(pvalue) >5), size = 4, color = 'black') + theme_light()  + xlim(-2,2)

only_aNSC <- medias[which(medias$type == 'aNSC Not DE' | medias$type == 'aNSC Downregulated'),]
p <- ggplot(only_aNSC, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel( size = 4, color = 'black') + theme_light()  + xlim(-2,2) + scale_fill_manual(values = dark_col[c(2,3)])

only_qNSC <- medias[which(medias$type == 'qNSC Not DE' | medias$type == 'qNSC Upregulated'),]
p <- ggplot(only_qNSC, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel( size = 4, color = 'black') + theme_light()  + xlim(-2,2) +scale_fill_manual(values = dark_col[c(3,1)])

only_DE <- medias[which(medias$type == 'qNSC Upregulated' | medias$type == 'aNSC Downregulated'),]
p <- ggplot(only_DE, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p2<- p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel( data = subset(only_DE, -log10(pvalue) > 2),size = 2.5, color = 'black') + theme_light()  + xlim(-2,2) +scale_fill_manual(values = dark_col[c(3,1)])
p2 + geom_hline(yintercept = 2, linetype = 'dashed', color = 'gray', size = 0.75, alpha = 1)


p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + theme_light()  + xlim(-2,2) +scale_fill_manual(values = dark_col[c(3,1)]) +  geom_hline(yintercept = 2, linetype = 'dashed', color = 'gray', size = 0.75)

bad_data <- medias[which(medias$type == 'qNSC Not DE' | medias$type == 'aNSC Not DE'),]
write.table(bad_data,'KEGG_not_DE_coculture.txt', sep = '\t')
##################Same with david
dark_col <- brewer.pal(3, name = 'Dark2')

all_bub_rw <- bind_rows(a_david_r, a_david_w, q_david_r, q_david_w)
all_bub_rw$type <- c(rep('aNSC Downregulated',nrow(a_david_r)), rep('aNSC Not DE', nrow(a_david_w)),
                     rep('qNSC Upregulated', nrow(q_david_r)), rep('qNSC Not DE', nrow(q_david_w)))

all_bub_rw %>% group_by(term,type) %>% summarise(means = mean(logFC, na.rm =T), pvalue = mean(adj_pval, na.rm = T), counter = n()) %>% as.data.frame()-> medias
p <- ggplot(medias, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel(data=subset(medias, -log10(pvalue) >5), size = 4, color = 'black') + theme_light()  + xlim(-2,2)

only_aNSC <- medias[which(medias$type == 'aNSC Not DE' | medias$type == 'aNSC Downregulated'),]
p <- ggplot(only_aNSC, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel( size = 4, color = 'black') + theme_light()  + xlim(-2,2) + scale_fill_manual(values = dark_col[c(2,3)])

only_qNSC <- medias[which(medias$type == 'qNSC Not DE' | medias$type == 'qNSC Upregulated'),]
p <- ggplot(only_qNSC, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel( size = 4, color = 'black') + theme_light()  + xlim(-2,2) +scale_fill_manual(values = dark_col[c(3,1)])

only_DE <- medias[which(medias$type == 'qNSC Upregulated' | medias$type == 'aNSC Downregulated'),]
p <- ggplot(only_DE, aes(x = means, y = -log10(pvalue), size = counter, fill = type, label = term))
p2<- p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel( data = subset(only_DE, -log10(pvalue) > 2),size = 2.5, color = 'black') + theme_light()  + xlim(-2,2) +scale_fill_manual(values = dark_col[c(3,1)])
p2 + geom_hline(yintercept = 2, linetype = 'dashed', color = 'gray', size = 0.75, alpha = 1)


p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + theme_light()  + xlim(-2,2) +scale_fill_manual(values = dark_col[c(3,1)]) +  geom_hline(yintercept = 3, linetype = 'dashed', color = 'darkred', size = 1.2)

bad_data <- medias[which(medias$type == 'qNSC Not DE' | medias$type == 'aNSC Not DE'),]
write.table(bad_data,'KEGG_not_DE_coculture.txt', sep = '\t')
#############################genes not in codega, what do they know?#########################

not_codega <- setdiff(row.names(cc_sig_lfc), row.names(m_codega))
not_codega_entrez <- bitr(not_codega, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

not_codega_kegg <- enrichKEGG(not_codega_entrez$ENTREZID, organism = 'mmu')
enrichMap(not_codega_kegg)
dotplot(not_codega_kegg, showCategory = 20)

test <- i_plot_bubbles(enrichment_obj = not_codega_kegg, results_obj = cc_sig_lfc, ID = 5, labels = T, title = 'not in codega')
test %>% group_by(term) %>% summarise(means = mean(logFC), pvalue = mean(adj_pval), counter = n()) %>% as.data.frame()-> medias

p <- ggplot(medias, aes(x = means, y = -log10(pvalue), size = counter, label = term, fill = means))
p + geom_point(alpha = 0.5, shape = 21, color = '#000000') + geom_text_repel( size = 4, color = 'black') + theme_light()  + xlim(-2,2)

#############################################p53 overlap###########################
p_53_targets <- read.csv('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/06_Analysis/02_Data/Supplemental_Table2.csv', header = T, stringsAsFactors = F)
venn.diagram(list(p_53_targets[,1],row.names(up_reg)), filename = 'up_cc_p53.tiff', category.names = c('P53 Target','upregulated'),hyper.test = T)

venn.diagram(list(p_53_targets[,1],row.names(down_reg)), filename = 'down_cc_p53.tiff', category.names = c('P53 Target','Downregulated'), hyper.test = T)

p_53_targets[,2] <- 'p53_target'
row.names(p_53_targets) <- p_53_targets[,1]
p_53_targets <- select(p_53_targets, 2)

pheatmap(select(as.data.frame(cc_sig_lfc), log2FoldChange), cluster_cols = F, annotation_row = p_53_targets ,show_rownames = F, clustering_method = 'ward.D2', color = col_pretty(12), breaks = seq(-6,6, by = 1))

inter <- intersect(row.names(p_53_targets), row.names(cc_sig_lfc))
only_de <- setdiff(row.names(cc_sig_lfc), row.names(p_53_targets))
only_p53 <- setdiff(row.names(p_53_targets), row.names(cc_sig_lfc))
no <- 17433-length(inter) - length(only_de) - length(only_p53)

test <- as.table(rbind(c(length(inter), length(only_de)), c(length(only_p53), no)))
fisher.test(test, alternative = 'less')
c
venn.diagram(list(row.names(p_53_targets),row.names(cc_sig_lfc)), filename = 'test.tiff', category.names = c('P53 Target','Regulated genes'), hyper.test = T, total.population = nrow(res_cc))


#############################################Use weird new GO terms##########
library('data.table')
library('reshape2')
go_codega <- read.csv('../codega_GO_summary.csv', header =T, strings = F)
#NEED to remove repeated entries
unique_codega <- unique(setDT(go_codega)) #This is fucking magic

log2_fsc <- rownames_to_column(log2_fsc)
log2_fsc$mean_lfc <- rowMeans(log2_fsc[,2:4], na.rm = T)


go_lfc <- left_join(unique_codega,log2_fsc, by = c('Gene' = 'rowname'))

go_lfc$DE <- NA
go_lfc[which(go_lfc$Gene %in% row.names(cc_genes)),]$DE <- 'DE'
go_lfc[is.na(go_lfc$DE),]$DE <- 'Not DE'


go_lfc %>% group_by(GO, Cell.Type) %>% summarise(means1 = mean(mean_lfc, na.rm = T), broad = unique(Broad), counter = sum(DE == 'DE')/n()) %>% as.data.frame()-> avg_GO

head(avg_GO)
avg_GO$percen <- avg_GO$means1*avg_GO$counter
test <- melt(avg_GO, id.vars = c('GO','broad','Cell.Type','counter'))

p <- ggplot(test, aes( x = GO, y =value, fill = variable, group = broad))
p + geom_bar(stat = 'identity',position=position_dodge()) + theme_light() + theme(axis.text.x =element_text(size  = 10,angle = 45, hjust = 1,vjust = 1))  + scale_color_distiller(palette = 'BuGn', direction = 1)+coord_flip() + facet_wrap(~Cell.Type, scales = 'free_y', ncol = 1, strip.position = 'left') + ylim(-1.5,1.5)


###########################pathview
res_cc_entrez <- bitr(row.names(cc_sig_lfc), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
res_entrez_fsc <- merge(res_cc_entrez, as.data.frame(cc_sig_lfc), by.x = 'SYMBOL', by.y = 'row.names')

test <- res_entrez_fsc$log2FoldChange
names(test) <- res_entrez_fsc$ENTREZID

t<- pathview(gene.data = test, pathway.id = 'mmu03320', species = 'mmu', limit = 5)
t<- pathview(gene.data = test, pathway.id = 'mmu00071', species = 'mmu', limit = 2)


pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],species = "hsa", out.suffix = "gse16873", kegg.native = T)

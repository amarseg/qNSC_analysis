################functions for Simona's data processing
i_plot_bubbles <- function(enrichment_obj, results_obj, labels, ID, title)
{
	require('GOplot')
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

join_bubble_objects <- function(bub1, bub2)
{
	all_bub <- bind_rows(bub1, bub2) #Average log change for plotting
	
	all_bub %>% group_by(term) %>% summarise(means = mean(logFC), pvalue = mean(adj_pval), counter = n()) %>% as.data.frame()-> medias
	return(medias)
	
}

log2_fold_change <- function(exprs_matrix, reference, samples)
{
	exprs_matrix <- exprs_matrix + 0.0001
	ref_mat <- exprs_matrix[,reference]
	samp_mat <- exprs_matrix[,samples]
	
	fc <- samp_mat/ref_mat
	log2fc <- log2(fc)
	return(log2fc)
	
}

comparison_and_enrichment <- function(de1,de2, enrichfun = 'enrichKEGG', k = 8, labels, plot_title)
{
	require(dplyr)
	require('clusterProfiler')
	pdf(plot_title)
	col_pretty <- colorRampPalette(c("blue", "azure1", "red"))
	#merge the two objects
	all_de <- merge(as.data.frame(de1), as.data.frame(de2), by = 'row.names')
	only_log2 <- dplyr::select(all_de, log2FoldChange.x, log2FoldChange.y)
	row.names(only_log2) <- all_de$Row.names
	colnames(only_log2) <- labels
	heatmap_obj <- pheatmap(only_log2, kmeans_k = k, breaks = seq(-5,5,by = 1), color = col_pretty(10))
	dev.off()
	cl <- heatmap_obj$kmeans$cluster
	cl_entrez <- bitr(names(cl), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
	cl_profiler <- merge(cl_entrez, as.data.frame(cl), by.x = 'ALIAS', by.y = 'row.names')
	compare <- compareCluster(data = cl_profiler, ENTREZID ~ cl, fun = enrichfun, org = 'mmu' )
	t <- dotplot(compare)
	ggsave(t, filename = paste0('enrichment_',plot_title))
	
}

comparison_GO_terms <- function(bub1, bub2, de1, de2)
{

	
	de1 <- rownames_to_column(as.data.frame(de1))
	de1_ids <- bitr(de1$rowname, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
	de1_obj <- inner_join(de1, de1_ids, by = c('rowname' = 'ALIAS'))
	
	de2 <- rownames_to_column(as.data.frame(de2))
	de2_ids <- bitr(de1$rowname, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
	de2_obj <- inner_join(de2, de2_ids, by = c('rowname' = 'ALIAS'))
	
	bub1$sample <- 'sample1'
	bub2$sample <- 'sample2'
	
	test <- rbind(bub1,bub2)
	test2 <- left_join(test, de1_obj, by = c('genes' = 'ENTREZID'))
	test3 <- left_join(test2, de2_obj, by = c('genes' = 'ENTREZID'))
	
	test3 %>% group_by(term) %>% summarise(means1 = mean(log2FoldChange.x), means2 = mean(log2FoldChange.y), pval = mean(adj_pval)) %>% as.data.frame()-> medias1
	
	medias1$dataset <- NA
	medias1[which(medias1$term %in% bub1$term & medias1$term %in% bub2$term),]$dataset <- 'Both sets'
	medias1[which(medias1$term %in% bub1$term & !(medias1$term %in% bub2$term)),]$dataset <- 'Set 1'
	medias1[which(!(medias1$term %in% bub1$term) & medias1$term %in% bub2$term),]$dataset <- 'Set 2'
	
	return(medias1)
}


category_plotter <- function(category_name)
{
	require('tidyverse')
	all_enrich <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/enrichments_table.txt', header = T, strings = F)
	all_enrich$plot.pval <- -log10(all_enrich$p.adjust)
	cat_enrichment <- all_enrich[base::grep(all_enrich$Description, pattern = category_name, ignore.case = T),]
	p <- ggplot(subset(cat_enrichment, Comparison %in% c('coculture_vs_control','p53_KOcoculture_vs_flox_coculture')), aes(Comparison,Description ))
	p2<- p + geom_tile(aes(fill = plot.pval)) + facet_grid(database~direction, drop = T) + theme_light() + theme(axis.text.x =element_text(size  = 10,angle = 45, hjust = 1,vjust = 1))
	
	
	print(p2)
	return(p2)
}

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
	
	all_enrich <- rbind(up_GO, down_GO, up_david, down_david, down_KEGG, up_KEGG)	
	all_enrich$Comparison <- comparison
	return(all_enrich)
}

gene_extractor <- function(GO_term, comparison, direction)
{
	require('clusterProfiler')
	require('Org.Mm.db.eg')
	enrichment_data <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/Simona_project/enrichments_table.txt', header = T, strings = F)
	select_data <- enrichment_data[which(enrichment_data$ID == GO_term & enrichment_data$Comparison == comparison & enrichment_data$direction == direction	),]
	
	ids <- as.vector(strsplit(select_data$geneID, split = '/', fixed = T))
	symbol <- bitr(ids, fromType = 'ENTREZID', toType = 'SYMBOL', org.db = org.Mm.eg.db)
	return(symbol)
}

heatmapper <- function(gene_list)
{
  require('pheatmap')
  require('RColorBrewer')
  
  data <- read.delim('P:/Simona_project/supp1.txt', header = T, strings = F)
  
  log2_data <- data[,c(1,2,5,11,17)]
  row.names(log2_data) <- log2_data$Gene_names
  log2_data <- log2_data[,-1]
  colnames(log2_data) <- c('Coculture','p53cc_vs_flox_cc','p53KO_vs_wt_alone','p53KOcc_vs_p53KO_alone')
  
  pval_data <- data[,c(1,4,7,13,19)]
  row.names(pval_data) <- pval_data$Gene_names
  pval_data <- pval_data[,-1]
  
  de_mat <- pval_data  
  de_mat[de_mat < 0.05] <- 'DE' 
  de_mat[de_mat != 'DE'] <- 'Not_DE'
  colnames(de_mat) <- c('Coculture','p53cc_vs_flox_cc','p53KO_vs_control','p53KOcc_vs_p53KO_alone')
  
  
  subset_log2 <- log2_data[which(row.names(log2_data) %in% gene_list),]
  
  col_pretty <- colorRampPalette(c("blue", "azure1", "red"))
  
  ann_colors = list()
  ann_colors[[1]] <- c(DE = 'black', Not_DE = 'white')
  for(i in 2:4){ann_colors[[i]] <- ann_colors[[1]]}
  names(ann_colors) <- colnames(de_mat)
  
  breaks <- seq(-ceiling(max(subset_log2)+0.5) , ceiling(max(subset_log2)+0.5),1)
  n = length(breaks)
  pheatmap(subset_log2[,c(1,3,2,4)], cluster_cols = F, annotation_row = de_mat[,c(4,2,3,1)], display_numbers = T, color = col_pretty(n), breaks = breaks, annotation_colors = ann_colors, show_colnames = T, annotation_names_row = T)

}
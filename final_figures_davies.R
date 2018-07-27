#This is a file that hopes to aggregate all the code to make figures for Tim's paper
################Input and libraries###########
#
library('tidyverse')
library('clusterProfiler')
library('DESeq2')
library('pheatmap')
library('IHW')
library('org.Mm.eg.db')
library('RDAVIDWebService')
library('RColorBrewer')
library('data.table')
library('reshape2')
library('wesanderson')
library('lettercase')
source('../figure_2/useful_funs.R')
library('ggrepel')
setwd('P:/Simona_project/final_figures/')

all_data <- read.delim('../supp1.txt')
all_enrich <- read.delim('../enrichments_table.txt')
col_pretty <- colorRampPalette(c("blue", "azure1", "red"))
################Boadilla and Codega##############
codega_data <- read.csv(file = 'P:/Simona_project/06_Analysis/02_Data/aNSC_qNSC_profiles.csv', header = T)

codega_data[codega_data == ""] <- NA
m_codega <- gather(codega_data,na.rm = T)
#I remove duplicates, how can this people upload a gene set with duplicates?!
names_split <- strsplit(m_codega$value, split = ' /// ', fixed = T)
good_names <- unlist(lapply(names_split, function(x){x[1]}))
m_codega$value <- good_names
m_codega <- m_codega[!duplicated(m_codega$value),]

row.names(m_codega) <- m_codega$value
m_codega <- dplyr::select(m_codega, key)

boadilla_data <- read.delim(file = 'P:/Simona_project/Boadilla_supp/bobadilla_clusters.txt', header = T, stringsAsFactors = F)
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


#################Heatmap codega y boadilla#############
seq_data <- load('P:/Simona_project/06_Analysis/02_Data/DATASET_070717.rda')
seq_data <- DATR
row.names(seq_data) <- seq_data$ID
seq_data <- seq_data[,-1]
#sample details
pheno_seq <- read.delim('P:/Simona_project/06_Analysis/02_Data/Samples_Davies_040717.txt', stringsAsFactors = F)
row.names(pheno_seq) <- pheno_seq$ID

p_adj_thr <- 0.05
res_cc <- all_data[which(all_data$padj_coculture_vs_control < p_adj_thr),]
fsc_cc <- (all_data[,c(23,27,31)]+0.0001)/(all_data[,c(22,26,30)]+0.0001)
row.names(fsc_cc) <- all_data$Gene_names
fsc_cc <- log2(fsc_cc + 0.00001)
annot_col <- list(CellType = c(qNSC = '#1B9E77', aNSC = "#D95F02"))

fsc_cc_DE <- fsc_cc[which(row.names(fsc_cc) %in% res_cc$Gene_names),]
pheatmap(fsc_cc_DE, show_rownames = F, color = col_pretty(30), annotation_row = m_all_id, 
         annotation_colors = annot_col[1],show_colnames = F, clustering_method = 'ward.D2', breaks = seq(-15,15,1))


#############################Barplot codega################
up_reg <- res_cc[which(res_cc$log2FoldChange_coculture_vs_control > 0),]
down_reg <- res_cc[which(res_cc$log2FoldChange_coculture_vs_control < 0),]
codega_up_down <- m_codega
codega_up_down$direction <- 'NA'
codega_up_down[which(row.names(codega_up_down) %in% up_reg$Gene_names & codega_up_down$key == 'qNSC'),]$direction <- 'Right Direction'
codega_up_down[which(row.names(codega_up_down) %in% down_reg$Gene_names & codega_up_down$key == 'aNSC'),]$direction <- 'Right Direction'

codega_up_down[which(row.names(codega_up_down) %in% down_reg$Gene_names & codega_up_down$key == 'qNSC'),]$direction <- 'Wrong Direction'
codega_up_down[which(row.names(codega_up_down) %in% up_reg$Gene_names & codega_up_down$key == 'aNSC'),]$direction <- 'Wrong Direction'

codega_up_down[which(codega_up_down$direction == 'NA'),]$direction <- 'Not Differentially regulated'

codega_up_down$direction <- factor(codega_up_down$direction, levels = c('Wrong Direction','Not Differentially regulated','Right Direction'))
write.table(codega_up_down, 'regulation_codega.txt', sep = '\t')


##########################Bubble plots KEGG#################


comparing_codega <- rownames_to_column(codega_up_down)

entrez_id <- bitr(comparing_codega$rowname, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

codega_entrez <- inner_join(comparing_codega, entrez_id, by = c('rowname' = 'ALIAS'))
right <- codega_entrez[codega_entrez$key == 'qNSC' & codega_entrez$direction == 'Right Direction',]
wrong <- codega_entrez[codega_entrez$key == 'qNSC' & codega_entrez$direction == 'Not Differentially regulated',]

r_w_list <- list(right$ENTREZID, wrong$ENTREZID)
names(r_w_list) <- c('right','Wrong')
m_kegg_comp <- compareCluster(r_w_list, 'enrichKEGG', organism = 'mmu')

right_a <- codega_entrez[codega_entrez$key == 'aNSC' & codega_entrez$direction == 'Right Direction',]
wrong_a <- codega_entrez[codega_entrez$key == 'aNSC' & codega_entrez$direction == 'Not Differentially regulated',]
r_w_list_a <- list(right_a$ENTREZID, wrong_a$ENTREZID)
names(r_w_list_a) <- c('right','Wrong')
m_kegg_comp_a <- compareCluster(r_w_list_a, 'enrichKEGG', organism = 'mmu')

kegg_df <- as.data.frame(m_kegg_comp)
kegg_df_a <- as.data.frame(m_kegg_comp_a)

res_cc <- res_cc[,c(1,2,3,4)]
row.names(res_cc) <- res_cc$Gene_names
res_cc <- res_cc[,-1]

kegg_df$Cluster <- as.factor(kegg_df$Cluster)
qNSC_kegg <- base::split(kegg_df, f = kegg_df$Cluster)
right_q_kegg <-i_plot_bubbles(qNSC_kegg[[1]], res_cc, 1, F, title = 'qNSC_right')
wrong_q_kegg <- i_plot_bubbles(qNSC_kegg[[2]], res_cc, 1, F, title = 'qNSC_wrong')

aNSC_kegg <- base::split(kegg_df_a, f = kegg_df_a$Cluster)
right_a_kegg <- i_plot_bubbles(aNSC_kegg[[1]],res_cc , 1, F, title = 'aNSC_righ')
wrong_a_kegg <- i_plot_bubbles(aNSC_kegg[[2]], res_cc, 1, F)

dark_col <- brewer.pal(3, name = 'Dark2')

all_bub_rw <- bind_rows(wrong_a_kegg, right_a_kegg, right_q_kegg, wrong_q_kegg)
all_bub_rw$type <- c(rep('aNSC Not DE',nrow(wrong_a_kegg)), rep('aNSC Downregulated', nrow(right_a_kegg)),
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

##########################IPA figures Upstream Coculture##############
setwd('P:/Simona_project/IPA_files/')
summarise_exp <- function( id_list, expr)
{
  id_list <- str_lowercase(id_list)
  id_list <- str_ucfirst(id_list)
  todo <- expr[which(expr[,1] %in% id_list),]
  x <- mean(todo[,2])
  return(x)
}
data <- all_data[,c(1,2)]
all_tfs <- read.delim('../IPA_files/TF_files/coculture_cutoff_data.txt', header = T, strings = F, skip = 2)
ids <- strsplit(all_tfs$Target.molecules.in.dataset, split = ',', fixed = T)
all_tfs$overlap <- unlist(lapply(ids, length))
all_tfs$mean_xp <- unlist(lapply(ids, summarise_exp,expr =  data))



all_codega <- read.delim('../IPA_files/TF_files/all_codega_tfs.txt', header = T, strings = F, skip = 2)
all_codega_all_tfs <- inner_join(all_tfs, all_codega, by = 'Upstream.Regulator')
all_codega_all_tfs[all_codega_all_tfs$Predicted.Activation.State.x == ' ',]$Predicted.Activation.State.x <- 'Inconclusive'
all_codega_all_tfs[all_codega_all_tfs$Predicted.Activation.State.y == ' ',]$Predicted.Activation.State.y <- 'Inconclusive'
all_codega_all_tfs<-all_codega_all_tfs[which(all_codega_all_tfs$Predicted.Activation.State.x != 'Inconclusive'),]

p <- ggplot(all_codega_all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap.x), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State.y))
p  + geom_point(alpha = 0.85, colour="black",pch=21, data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 4)) + geom_text_repel(data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 4),size = 2.5) + theme_light() + xlim(c(-5,5)) + facet_wrap(~Predicted.Activation.State.x, scale = 'free') + scale_fill_manual(values = wes_palette('BottleRocket2')[c(2,1,3)], name = 'Activation state in Codega')


p <- ggplot(all_codega_all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap.x), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State.y))
p  + geom_point(alpha = 0.85, colour="black",pch=21, data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2 & Predicted.Activation.State.y != 'Inconclusive')) + geom_text_repel(data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2),size = 2.5) + theme_light() + xlim(c(-5,5)) + facet_wrap(~Predicted.Activation.State.x, scale = 'free') + scale_fill_manual(values = wes_palette('BottleRocket2')[c(2,1,3)], name = 'Activation state in Codega')

################IPA figure KO upstream regulators

ko_cutoff_tfs <- read.delim('../IPA_files/TF_files/KO_cutoff_good.txt', header = T, strings = F, skip = 2)
#cutoff_tfs <- cutoff_tfs[which(cutoff_tfs$Predicted.Activation.State != " "),]
ids <- strsplit(ko_cutoff_tfs$Target.molecules.in.dataset, split = ',', fixed = T)
ko_cutoff_tfs$overlap <- unlist(lapply(ids, length))
ko_cutoff_tfs$mean_xp <- unlist(lapply(ids, summarise_exp,expr =  data))

all_codega_all_tfs <- inner_join(ko_cutoff_tfs, all_codega, by = 'Upstream.Regulator')
all_codega_all_tfs[all_codega_all_tfs$Predicted.Activation.State.x == ' ',]$Predicted.Activation.State.x <- 'Inconclusive'
all_codega_all_tfs[all_codega_all_tfs$Predicted.Activation.State.y == ' ',]$Predicted.Activation.State.y <- 'Inconclusive'

all_codega_all_tfs<-all_codega_all_tfs[which(all_codega_all_tfs$Predicted.Activation.State.x != 'Inconclusive'),]

p <- ggplot(all_codega_all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap.x), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State.y))
p  + geom_point(alpha = 0.85, colour="black",pch=21, data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2)) + geom_text_repel(data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2),size = 2.5) + theme_light() + xlim(c(-5,5)) + facet_wrap(~Predicted.Activation.State.x, scale = 'free') + scale_fill_manual(values = wes_palette('BottleRocket2')[c(2,1,3)], name = 'Activation state in Codega')

p <- ggplot(subset(all_codega_all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap.x), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State.y))
            p  + geom_point(alpha = 0.85, colour="black",pch=21, data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2 & Predicted.Activation.State.y != 'Inconclusive')) + geom_text_repel(data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2),size = 2.5) + theme_light() + xlim(c(-5,5)) + facet_wrap(~Predicted.Activation.State.x, scale = 'free') + scale_fill_manual(values = wes_palette('BottleRocket2')[c(2,1,3)], name = 'Activation state in Codega')
            
#####################Fatty acid heatmap
          
data <- read.delim('P:/Simona_project/supp1.txt', header = T, strings = F)
lipid_go <- read.delim('go_fatty_Acid.txt', header = F, strings = F)
lipid_go <- unique(lipid_go[,2])
lipid_go <- append(lipid_go,'Acad11')
            
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
colnames(de_mat) <- c('Coculture','p53cc_vs_flox_cc','p53KO_vs_wt_alone','p53KOcc_vs_p53KO_alone')
            
            
subset_log2 <- log2_data[which(row.names(log2_data) %in% lipid_go),]
          
col_pretty <- colorRampPalette(c("blue", "azure1", "red"))
            
ann_colors = list()
ann_colors[[1]] <- c(DE = 'black', Not_DE = 'white')
for(i in 2:4){ann_colors[[i]] <- ann_colors[[1]]}
names(ann_colors) <- colnames(de_mat)
            
breaks <- seq(-4.5 ,4.5,1)
n = length(breaks)
pheatmap(subset_log2[,c(1,4,2,3)], cluster_cols = F, annotation_row = de_mat[,c(3,2,4,1)], display_numbers = T, color = col_pretty(n), breaks = breaks, annotation_colors = ann_colors, show_colnames = T, annotation_names_row = T)
pheatmap(subset_log2[,c(1,4,2,3)], cluster_cols = F, annotation_row = de_mat[,c(3,2,4,1)], display_numbers = T, color = col_pretty(n), breaks = breaks, annotation_colors = ann_colors, show_colnames = F, annotation_names_row = F)


############################TF heatmaps

aNSC <- read.delim('../IPA_files/TF_files/aNSC_fold_change.txt', header = T, strings = F, skip = 2)
aNSC <- aNSC[,c(1,7)]
aNSC <- aNSC[aNSC$p.value.of.overlap < 0.05,]
qNSC <- read.delim('../IPA_files/TF_files/qNSC_fold_change.txt', header = T, strings = F, skip = 2)
qNSC <- qNSC[,c(1,7)]
qNSC <- qNSC[qNSC$p.value.of.overlap < 0.05,]
up <- read.delim('../IPA_files/TF_files/up_regulated_tfs.txt', header = T, strings = F, skip = 2)
up <- up[,c(1,7)]
up <- up[up$p.value.of.overlap < 0.05,]
down <- read.delim('../IPA_files/TF_files/down_regulated_tfs.txt', header = T, strings = F, skip = 2)
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


####################################Barplot with codega GO terms
all_data <- read.delim('../supp1.txt', header = T, strings = F)

coculture_data <- all_data[,1:4]
coculture_DE <- coculture_data[which(coculture_data$padj_coculture_vs_control < 0.05),]
KO_coculture_data <- all_data[,c(1,5:7)]
KO_coculture_DE <- KO_coculture_data[which(KO_coculture_data$padj_p53KO_coculture_vs_flox_coculture < 0.05),]


go_codega <- read.csv('../codega_GO_summary.csv', header =T, strings = F)
#NEED to remove repeated entries
unique_codega <- unique(setDT(go_codega)) #This is fucking magic

go_lfc <- left_join(unique_codega,coculture_data, by = c('Gene' = 'Gene_names'))
go_lfc <- left_join(go_lfc,KO_coculture_data, by = c('Gene' = 'Gene_names'))



go_lfc$DE_cc <- NA
go_lfc[which(go_lfc$Gene %in% coculture_DE$Gene_names),]$DE_cc <- 'DE'
go_lfc[is.na(go_lfc$DE_cc),]$DE_cc <- 'Not DE'

go_lfc$DE_KO <- NA
go_lfc[which(go_lfc$Gene %in% KO_coculture_DE$Gene_names),]$DE_KO <- 'DE'
go_lfc[is.na(go_lfc$DE_KO),]$DE_KO <- 'Not DE'

go_lfc %>% group_by(GO, Cell.Type) %>% summarise(mean_log2_fold_change_coculture = mean(log2FoldChange_coculture_vs_control, na.rm = T), broad = unique(Broad), counter_CC = sum(DE_cc == 'DE')/n(), mean_log2_fold_change_p53_KO_coculture = mean(log2FoldChange_p53KO_coculture_vs_flox_coculture, na.rm = T),counter_KO = sum(DE_KO == 'DE')/n() ) %>% as.data.frame()-> avg_GO

head(avg_GO)
avg_GO$proportion_DE_coculture <- avg_GO$mean_log2_fold_change_coculture*avg_GO$counter_CC
avg_GO$proportion_DE_KO_coculture <- avg_GO$mean_log2_fold_change_p53_KO_coculture*avg_GO$counter_KO

test <- melt(avg_GO, id.vars = c('GO','Cell.Type'), measure.vars = c('mean_log2_fold_change_coculture','mean_log2_fold_change_p53_KO_coculture','proportion_DE_coculture','proportion_DE_KO_coculture'))

test$comp <- NA
test[which(test$variable == 'mean_log2_fold_change_coculture' | test$variable == 'proportion_DE_coculture'),]$comp <- 'Coculture'
test[is.na(test$comp),]$comp <- 'KO_Coculture'

col <- wes_palette(n=4, name="Zissou1")
col_g <- col[c(2,1)]

p <- ggplot(test, aes( x = GO, y =value, fill = variable, group = comp))
p + geom_bar(stat = 'identity',position=position_dodge(), data = subset(test, comp == 'Coculture')) + theme_light() + theme(axis.text.x =element_text(size  = 10,angle = 45, hjust = 1,vjust = 1))  + scale_color_distiller(palette = 'BuGn', direction = 1)+coord_flip() + facet_wrap(~Cell.Type, scales = 'free_y', ncol = 1, strip.position = 'left') + ylim(-1.5,1.5)+ scale_fill_manual(values = col_g) 

col_g <- col[c(2,3,1,4)]

p <- ggplot(test, aes( x = GO, y =value, fill = variable, group = comp))
p + geom_bar(stat = 'identity',position=position_dodge()) + theme_light() + theme(axis.text.x =element_text(size  = 10,angle = 45, hjust = 1,vjust = 1))  + scale_color_distiller(palette = 'BuGn', direction = 1)+coord_flip() + facet_wrap(~Cell.Type, scales = 'free_y', ncol = 1, strip.position = 'left') + ylim(-1.5,1.5)+ scale_fill_manual(values = col_g) 

            
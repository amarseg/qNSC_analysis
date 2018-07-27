require('pheatmap')
require('RColorBrewer')
setwd('P:/Simona_project/')
library(extrafont)
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.22/bin/gswin64c.exe")


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

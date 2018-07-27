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
library('data.table')
library('reshape2')
library('wesanderson')
setwd('P:/Simona_project/figure_1')

###############Load data
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

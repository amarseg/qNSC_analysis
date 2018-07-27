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
library('gage')
library('GSVA')

setwd('P:/Simona_project/')
source('figure_2/useful_funs.R')

gmt_file <- readList('msigdb.v6.1.symbols.gmt')

supp <- read.delim('supp1.txt', header = T, strings = F)
cc_data <- supp[,c(1:4)]

cc_data$Gene_names <- toupper(cc_data$Gene_names)


data <- data.frame(cc_data$stat_coculture_vs_control, row.names = cc_data$Gene_names)

gsea_set <- gage(data, gset = gmt_file, ref = NULL, samp = NULL,rank.test=T)

great_gsea <- as.data.frame(gsea_set$greater)
gsea_sig <- great_gsea[which(great_gsea$q.val < 0.05),]

write.table(gsea_sig, 'gsea_overrepresented_coculture.txt', sep = '\t')

less_gsea <- as.data.frame(gsea_set$less)
less_sig <- less_gsea[which(less_gsea$q.val < 0.05),]

write.table(less_sig, 'gsea_underrepresented_coculture.txt', sep = '\t')

############compare with codega
codega <- read.delim('gsea_summart.txt', header = T,strings = F)

codega_in_vitro <- codega[which(codega$Description %in% row.names(gsea_sig)),]
write.table(codega_in_vitro,'overrepresented_in_common.txt', sep = '\t')

codega_in_vitro_less <- codega[which(codega$Description %in% row.names(less_sig)),]
write.table(codega_in_vitro_less,'underrepresented_in_common.txt', sep = '\t')

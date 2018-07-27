###############
rm(list = ls())
library('tidyverse')
library(clusterProfiler)
source('P:/Simona_project/figure_2/useful_funs.R')
library('org.Mm.eg.db')
library('ggrepel')

setwd('P:/Simona_project/')

supp <- read.delim('supp1.txt', header = T, strings = F)

data <- supp[,c(1,5:7)]
row.names(data) <- data$Gene_names

codega_data <- read.csv(file = 'P:/Simona_project/06_Analysis/02_Data/aNSC_qNSC_profiles.csv', header = T)
codega_data[codega_data == ""] <- NA
m_codega <- gather(codega_data,na.rm = T)
#I remove duplicates, how can this people upload a gene set with duplicates?!
names_split <- strsplit(m_codega$value, split = ' /// ', fixed = T)
good_names <- unlist(lapply(names_split, function(x){x[1]}))
m_codega$value <- good_names
m_codega <- m_codega[!duplicated(m_codega$value),]


sig_data <- data[which(data$padj_p53KO_coculture_vs_flox_coculture < 0.05),]

codega_data <- inner_join(sig_data, m_codega, by=c('Gene_names' = 'value'))

data <- data[,-1]
qNSC_data <- codega_data[codega_data$key == 'qNSC' & codega_data$log2FoldChange_p53KO_coculture_vs_flox_coculture < 0,]
qNSC_entrez <- bitr(qNSC_data$Gene_names, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
aNSC_data <- codega_data[codega_data$key == 'aNSC' & codega_data$log2FoldChange_p53KO_coculture_vs_flox_coculture > 0,]
aNSC_entrez <- bitr(aNSC_data$Gene_names, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)


kegg_q <- enrichKEGG(qNSC_entrez$ENTREZID, organism = 'mmu')
circle_q <- i_plot_bubbles(kegg_q, data, labels = 0, ID = F, title = 'qNSC')
kegg_a <- enrichKEGG(aNSC_entrez$ENTREZID, organism = 'mmu')
circle_a <- i_plot_bubbles(kegg_a, data, labels = 0, ID = F, title = 'qNSC')

t <- join_bubble_objects(circle_a, circle_q)

t$cell <- 'aNSC'
t[which(t$term %in% circle_q$term),]$cell <- 'qNSC'

p <- ggplot(t, aes(x = means, y = -log10(pvalue), size = counter, label = term, fill = cell))
p2<- p + geom_point(shape = 21,  alpha = 0.75) 
p2 +  geom_text_repel(size = 3) + scale_fill_brewer(palette = 'Dark2') + theme_light()

#now DAVID
kegg_q <- enrichDAVID(qNSC_entrez$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk', annotation = 'GOTERM_BP_4')
circle_q <- i_plot_bubbles(kegg_q, data, labels = 1, ID = F, title = 'qNSC')
kegg_a <- enrichDAVID(aNSC_entrez$ENTREZID, david.user = 'a.martinez-segura13@imperial.ac.uk', annotation = 'GOTERM_BP_4')
circle_a <- i_plot_bubbles(kegg_a, data, labels = 1, ID = F, title = 'qNSC')

t <- join_bubble_objects(circle_a, circle_q)

t$cell <- 'aNSC'
t[which(t$term %in% circle_q$term),]$cell <- 'qNSC'

p <- ggplot(t[t$cell == 'qNSC',], aes(x = means, y = -log10(pvalue), size = counter, label = term, fill = cell))
p2<- p + geom_point(shape = 21,  alpha = 0.75) 
p2 +  geom_text_repel(size = 3) + scale_fill_brewer(palette = 'Dark2') + theme_light() + facet_wrap(~cell, scales = 'free') + xlim(c(-1,1))


lipid_go <- read.delim('../go_fatty_Acid.txt', header = F, strings = F)
lipid_go <- unique(lipid_go[,2])
lipid_go <- append(lipid_go,'Acad11')
heatmapper(lipid_go)

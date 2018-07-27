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
setwd('P:/Simona_project/')
source('P:/Simona_project/06_Analysis/02_Data/Grouping.R')
source('figure_2/useful_funs.R')

supp <- read.delim('supp1.txt', header = T, strings = F)
supp_cc <- supp[,c(1:4)]

supp_cc$Sample <- 'Coculture_WT'

p53_cc <- supp[,c(1,5:7)]  
p53_cc$Sample <- 'Coculture_p53'
colnames(p53_cc) <- c('Gene_Name', 'Log2','Stat','Pval','Sample')
colnames(supp_cc) <- c('Gene_Name', 'Log2','Stat','Pval','Sample')

test <- base::rbind(supp_cc, p53_cc)

library('data.table')
library('reshape2')
go_codega <- read.csv('codega_GO_summary.csv', header =T, strings = F)
#NEED to remove repeated entries
unique_codega <- unique(setDT(go_codega)) #This is fucking magic


go_lfc <- left_join(unique_codega,test, by = c('Gene' = 'Gene_Name'))

go_lfc$DE <- NA
go_lfc[which(go_lfc$Pval < 0.05),]$DE <- 'DE'
go_lfc[is.na(go_lfc$DE),]$DE <- 'Not DE'


go_lfc %>% group_by(GO, Cell.Type,Sample) %>% summarise(means1 = mean(Log2, na.rm = T), broad = unique(Broad), counter = sum(DE == 'DE')/n()) %>% as.data.frame()-> avg_GO

head(avg_GO)
avg_GO$percen <- avg_GO$means1*avg_GO$counter
test2 <- melt(avg_GO, id.vars = c('GO','broad','Cell.Type','counter','Sample'))

p <- ggplot(test2, aes( x = GO, y =value, group = Sample, fill))
p + geom_bar(stat = 'identity',position=position_dodge()) + theme_light() + theme(axis.text.x =element_text(size  = 10,angle = 45, hjust = 1,vjust = 1))  + scale_color_distiller(palette = 'BuGn', direction = 1)+coord_flip() + facet_wrap(~Cell.Type, scales = 'free_y', ncol = 1, strip.position = 'left') + ylim(-1.5,1.5)


###########extra jak/stat and butaonate pathways
jak <- read.delim('P:/Simona_project/jak_stat_kegg.txt', header = F)

jak_entrez <- bitr(jak$V1, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = org.Mm.eg.db)

heatmapper(jak_entrez$SYMBOL)

but <- read.delim('P:/Simona_project/butanoate_kegg.txt', header = F)
but_entrez <- bitr(but$V1, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = org.Mm.eg.db)

heatmapper(but_entrez$SYMBOL)

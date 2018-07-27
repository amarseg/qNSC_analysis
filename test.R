library('DESeq2')

supp <- read.delim('P:/Simona_project/supp1.txt', header = T, strings = F)

load('P:/Simona_project/figure_2/grouped_DESeq.RData')

cc <- results(grouped_ds, contrast = c('group','floxendo','floxcontrol_endo'))

cc <- as.data.frame(cc)

test <- supp[,c(1:4)]

test$Gene_names == row.names(cc)

plot(test$log2FoldChange_coculture_vs_control, cc$log2FoldChange)
plot(test$padj_coculture_vs_control, cc$padj)

vst_ds <- vst(grouped_ds)
plotPCA(vst_ds, intgroup = 'group')

sig <- test[which(test$padj_coculture_vs_control <0.05),]
sig2 <- cc[which(cc$padj <0.05),]

#########IPA interface for plotting
rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library('ggrepel')
library('lettercase')
library('wesanderson')

setwd('P:/Simona_project/IPA_files/')
summarise_exp <- function( id_list, expr)
{
  id_list <- str_lowercase(id_list)
  id_list <- str_ucfirst(id_list)
  todo <- expr[which(expr[,1] %in% id_list),]
  x <- mean(todo[,2])
  return(x)
}

supp <- read.delim('../supp1.txt', strings = F, header =T)
data <- supp[,c(1,2)]

tfs_up <- read.delim('TF_files/up_regulated_tfs.txt', header = T, strings = F, skip = 2)
tfs_up <- tfs_up[,-6]
ids <- strsplit(tfs_up$Target.molecules.in.dataset, split = ',', fixed = T)
tfs_up$overlap <- unlist(lapply(ids, length))
tfs_up$mean_xp <- unlist(lapply(ids, summarise_exp,expr =  data))
p <- ggplot(tfs_up, aes(x = mean_xp, y = -log10(p.value.of.overlap), size = overlap, label = Upstream.Regulator, colour = Predicted.Activation.State))
p  + geom_text(data = subset(tfs_up, -log10(p.value.of.overlap) > 3))

tfs_down <- read.delim('TF_files/down_regulated_tfs.txt', header = T, strings = F, skip = 2)

ids <- strsplit(tfs_down$Target.molecules.in.dataset, split = ',', fixed = T)
tfs_down$overlap <- unlist(lapply(ids, length))
tfs_down$mean_xp <- unlist(lapply(ids, summarise_exp,expr =  data))
p <- ggplot(tfs_down, aes(x = mean_xp, y = -log10(p.value.of.overlap), size = overlap, label = Upstream.Regulator, colour = Predicted.Activation.State))
p  + geom_text(data = subset(tfs_down, -log10(p.value.of.overlap) > 3))



all_tfs <- read.delim('TF_files/coculture_cutoff_data.txt', header = T, strings = F, skip = 2)
all_tfs <- all_tfs[which(all_tfs$Predicted.Activation.State != " "),]
ids <- strsplit(all_tfs$Target.molecules.in.dataset, split = ',', fixed = T)
all_tfs$overlap <- unlist(lapply(ids, length))
all_tfs$mean_xp <- unlist(lapply(ids, summarise_exp,expr =  data))
p <- ggplot(all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State, color = Predicted.Activation.State))
p  + geom_text(data = subset(all_tfs, -log10(p.value.of.overlap) > 5))

pal <- as.vector(wes_palette('Darjeeling1'))
pal <- pal[c(3,2,1)]
p  + geom_point(alpha = 0.75, colour="black",pch=21) + geom_text_repel(data = subset(all_tfs, -log10(p.value.of.overlap) > 5), size = 3) + theme_light() + xlim(c(-5,5)) + scale_fill_manual(values = pal) + scale_color_manual(values = pal)


cutoff_tfs <- read.delim('TF_files/all_tfs_0.01_cutoff.txt', header = T, strings = F, skip = 2)
cutoff_tfs <- cutoff_tfs[which(cutoff_tfs$Predicted.Activation.State != " "),]
ids <- strsplit(cutoff_tfs$Target.molecules.in.dataset, split = ',', fixed = T)
cutoff_tfs$overlap <- unlist(lapply(ids, length))
cutoff_tfs$mean_xp <- unlist(lapply(ids, summarise_exp,expr =  data))
p <- ggplot(cutoff_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap), size = overlap, label = Upstream.Regulator, colour = Predicted.Activation.State, fill = Predicted.Activation.State))
pal <- as.vector(wes_palette('Darjeeling1'))
pal <- pal[c(2,1)]
p  + geom_point(data = subset(cutoff_tfs, -log10(p.value.of.overlap) > 2),alpha = 0.75, colour="black",pch=21) + geom_text_repel(data = subset(cutoff_tfs, -log10(p.value.of.overlap) > 2), size = 3) + theme_light() + xlim(c(-3,3)) + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + facet_wrap(~Predicted.Activation.State)
p  + geom_point(,alpha = 0.75, colour="black",pch=21) + geom_text_repel(size = 3) + theme_light() + xlim(c(-3,3)) + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + facet_wrap(~Predicted.Activation.State)


ko_cutoff_tfs <- read.delim('TF_files/KO_cutoff_good.txt', header = T, strings = F, skip = 2)
#cutoff_tfs <- cutoff_tfs[which(cutoff_tfs$Predicted.Activation.State != " "),]
ids <- strsplit(ko_cutoff_tfs$Target.molecules.in.dataset, split = ',', fixed = T)
ko_cutoff_tfs$overlap <- unlist(lapply(ids, length))
ko_cutoff_tfs$mean_xp <- unlist(lapply(ids, summarise_exp,expr =  data))
p <- ggplot(ko_cutoff_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap), size = overlap, label = Upstream.Regulator, colour = Predicted.Activation.State, fill = Predicted.Activation.State))
pal <- as.vector(wes_palette('Darjeeling1'))
pal <- pal[c(3,2,1)]
p  + geom_point(data = subset(ko_cutoff_tfs, -log10(p.value.of.overlap) > 3),alpha = 0.75, colour="black",pch=21) + geom_text_repel(data = subset(ko_cutoff_tfs, -log10(p.value.of.overlap) > 3), size = 3) + theme_light() + xlim(c(-3,3)) + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + facet_wrap(~Predicted.Activation.State)
p  + geom_point(alpha = 0.75, colour="black",pch=21) + geom_text_repel(size = 3) + theme_light() + xlim(c(-3,3)) + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + facet_wrap(~Predicted.Activation.State)

############Same with pathways
path_up <- read.delim('TF_files/pathways_up_good.txt', header = T, skip = 2, strings =F)
ids <- strsplit(path_up$Molecules, split = ',', fixed = T)
path_up$overlap <- unlist(lapply(ids, length))
path_up$avg <- unlist(lapply(ids, summarise_exp, expr = data))

path_down <- read.delim('TF_files/pathways_up.txt', header = T, skip = 2, strings = F)
ids <- strsplit(path_down$Molecules, split = ',', fixed = T)
path_down$overlap <- unlist(lapply(ids, length))
path_down$avg <- unlist(lapply(ids, summarise_exp, expr = data))

pal <- wes_palette("Zissou1", 100, type = "continuous")
all_paths <- dplyr::bind_rows(path_up, path_down)
p <- ggplot(all_paths, aes(x = avg, y = X.log.p.value., size = overlap, label = Ingenuity.Canonical.Pathways, fill = avg))
p  + geom_point(alpha = 0.75, colour="black",pch=21) + geom_text_repel(data = subset(all_paths, X.log.p.value. > 7), size = 3) + theme_light() + xlim(c(-5,5)) + scale_fill_gradientn(colours = pal)

p  + geom_point(data = subset(all_paths, X.log.p.value. > 6),alpha = 0.75, colour="black",pch=21) + geom_text_repel(data = subset(all_paths, X.log.p.value. > 6.5), size = 3) + theme_light() + xlim(c(-5,5)) + scale_fill_gradientn(colours = pal)

#####Now I need to do the fucking comparison with the codega dataset 


all_tfs <- read.delim('TF_files/coculture_cutoff_data.txt', header = T, strings = F, skip = 2)
ids <- strsplit(all_tfs$Target.molecules.in.dataset, split = ',', fixed = T)
all_tfs$overlap <- unlist(lapply(ids, length))
all_tfs$mean_xp <- unlist(lapply(ids, summarise_exp,expr =  data))


all_codega <- read.delim('TF_files/all_codega_tfs.txt', header = T, strings = F, skip = 2)
all_codega_all_tfs <- inner_join(all_tfs, all_codega, by = 'Upstream.Regulator')
all_codega_all_tfs[all_codega_all_tfs$Predicted.Activation.State.x == ' ',]$Predicted.Activation.State.x <- 'Inconclusive'
all_codega_all_tfs[all_codega_all_tfs$Predicted.Activation.State.y == ' ',]$Predicted.Activation.State.y <- 'Inconclusive'
all_codega_all_tfs<-all_codega_all_tfs[which(all_codega_all_tfs$Predicted.Activation.State.x != 'Inconclusive'),]

p <- ggplot(all_codega_all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap.x), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State.y))
p  + geom_point(alpha = 0.85, colour="black",pch=21, data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 4)) + geom_text_repel(data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 4),size = 2.5) + theme_light() + xlim(c(-5,5)) + facet_wrap(~Predicted.Activation.State.x, scale = 'free') + scale_fill_manual(values = wes_palette('BottleRocket2')[c(2,1,3)], name = 'Activation state in Codega')


p <- ggplot(all_codega_all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap.x), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State.y))
p  + geom_point(alpha = 0.85, colour="black",pch=21, data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2 & Predicted.Activation.State.y != 'Inconclusive')) + geom_text_repel(data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2),size = 2.5) + theme_light() + xlim(c(-5,5)) + facet_wrap(~Predicted.Activation.State.x, scale = 'free') + scale_fill_manual(values = wes_palette('BottleRocket2')[c(2,1,3)], name = 'Activation state in Codega')
            


all_codega_all_tfs <- inner_join(ko_cutoff_tfs, all_codega, by = 'Upstream.Regulator')
all_codega_all_tfs[all_codega_all_tfs$Predicted.Activation.State.x == ' ',]$Predicted.Activation.State.x <- 'Inconclusive'
all_codega_all_tfs[all_codega_all_tfs$Predicted.Activation.State.y == ' ',]$Predicted.Activation.State.y <- 'Inconclusive'

all_codega_all_tfs<-all_codega_all_tfs[which(all_codega_all_tfs$Predicted.Activation.State.x != 'Inconclusive'),]

p <- ggplot(all_codega_all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap.x), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State.y))
p  + geom_point(alpha = 0.85, colour="black",pch=21, data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2)) + geom_text_repel(data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2),size = 2.5) + theme_light() + xlim(c(-5,5)) + facet_wrap(~Predicted.Activation.State.x, scale = 'free') + scale_fill_manual(values = wes_palette('BottleRocket2')[c(2,1,3)], name = 'Activation state in Codega')

p <- ggplot(subset(all_codega_all_tfs, aes(x = mean_xp, y = -log10(p.value.of.overlap.x), size = overlap, label = Upstream.Regulator, fill = Predicted.Activation.State.y))
p  + geom_point(alpha = 0.85, colour="black",pch=21, data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2 & Predicted.Activation.State.y != 'Inconclusive')) + geom_text_repel(data = subset(all_codega_all_tfs, -log10(p.value.of.overlap.x) > 2),size = 2.5) + theme_light() + xlim(c(-5,5)) + facet_wrap(~Predicted.Activation.State.x, scale = 'free') + scale_fill_manual(values = wes_palette('BottleRocket2')[c(2,1,3)], name = 'Activation state in Codega')

library(extrafont)
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.22/bin/gswin64c.exe")

embed_fonts('coculture_IPA_final.pdf', outfile = 'coculture_IPA_final_Emb.pdf')
embed_fonts('KO_IPA_plot.pdf', outfile = 'KO_IPA_final_Emb.pdf')


library("ggplot2")
library("reshape2")
library("aplot")
library("tidyr")

rm(list=ls())
pathway <- read.csv("KEGG_diff_d_compare_CK_for_heatmap.csv",header = T, sep = ",", check.names = F)

pathway_longdata <- melt(pathway)
pathway_order <- unique(as.character(pathway_longdata$pathway))
# draw heatmap
p_cor <- ggplot(pathway_longdata, aes(factor(pathway, levels = pathway_order), variable)) + 
  geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(name="p value", low = "red",high = "white") +
  theme(axis.text.x = element_text(hjust = 1, angle = 75))+
  coord_fixed(ratio=1)+
  labs(x = "Pathway", y = "Group", title = "p-value")+
  theme(plot.title = element_text(size = 13,hjust = 0.5))

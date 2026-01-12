
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationHub")
library("dplyr")
library(ggplot2)
library(ggrepel)
library("AnnotationHub")
library("clusterProfiler") 

term2gene <- read.csv("././term2gene_all.csv",header=T,sep = ",")

term2name <- read.csv("././term2name_all.csv",header = T,sep = ",")

gene <- read.csv("CK_vs_Y0_75_2d.deseq2_down_geneid.csv",header = F,sep = ",")
gene <- as.factor(gene$V1)

x <- enricher(gene,TERM2GENE = term2gene,TERM2NAME = term2name,pvalueCutoff = 0.05,
              pAdjustMethod = "none",qvalueCutoff = 1)
head(x, 10)

write.csv(x, "KEGG_CK_vs_Y0_75_2d.deseq2_down.csv", quote = F, row.names = F)
pvalue_ <- pvalue<0.05
dotplot(x, x="Count", color="pvalue")

p <- ggplot(KEGG_dataset,aes(x=GeneRatio,y=Description,colour=pvalue,size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("GeneRatio")+
  labs(color = expression(pvalue))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size = 14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))
p




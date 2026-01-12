# 加载必要的库
library(ggplot2)
BiocManager::install("ggalt")
library(factoextra)
library(FactoMineR)
library(ggbiplot)

# 读取数据
data_path <- "e:/Nostoc盐胁迫/RNA_seq_1week/salt-stress_R/PCA/transcript.tpm.matrix_all.csv"
data <- read.csv(data_path, header = TRUE, row.names = 1)

# 转置数据矩阵，使样本为行，基因为列
data_t <- t(data)

# 检查并处理缺失值或无穷值
data_t[is.na(data_t)] <- 0  # 将NA替换为0
data_t[!is.finite(data_t)] <- 0  # 将无穷值替换为0

# 移除全为零的行
data_t <- data_t[rowSums(data_t) > 0, ]

# 移除方差为零的列（常数列）
data_t <- data_t[, apply(data_t, 2, var) > 0]

# 数据标准化前先进行对数转换（处理基因表达数据常用方法）
# 添加一个小常数避免log(0)的问题
data_log <- log2(data_t + 1)

# 正确提取分组信息
# 根据样本名称的实际格式创建分组
sample_names <- rownames(data_t)
groups <- vector("character", length(sample_names))


for (i in 1:length(sample_names)) {
  name <- sample_names[i]
  if (grepl("^CK", name)) {
    groups[i] <- "CK"
  } else if (grepl("^Y0M_2d", name)) {
    groups[i] <- "Y0M_2d"
  } else if (grepl("^Y0_25", name)) {
    groups[i] <- "Y0_25"
  } else if (grepl("^Y0_5_1", name)) {
    groups[i] <- "Y0_5_1d"
  } else if (grepl("^Y0_5_2", name)) {
    groups[i] <- "Y0_5_2d"
  } else if (grepl("^Y0_5_4", name)) {
    groups[i] <- "Y0_5_4d"
  } else if (grepl("^Y0_75", name)) {
    groups[i] <- "Y0_75"
  } else {
    groups[i] <- name  # 如果没有匹配，保留原始名称
  }
}

# 确保分组是因子类型，并指定顺序
groups <- factor(groups, levels = c("CK", "Y0M_2d", "Y0_25",  "Y0_5_1d", "Y0_5_2d", "Y0_5_4d", "Y0_75"))

rld <- rlogData(data_t)
plotPCA(rld, intgroup=c("condition","period"),returnData =
          F)

# 使用相关系数矩阵而不是协方差矩阵进行PCA
pca_result1 <- prcomp(data_log, center = TRUE, scale. = TRUE)
write.table(pca_result1$rotation, file="PC.xls", quote=F, sep = "\t") #输出新表 
write.table(predict(pca_result), file="newTab1.xls", quote=F, sep = "\t") #输出PC比重 
pca.sum=summary(pca_result) 
write.table(pca.sum$importance, file="importance.xls", quote=F, sep = "\t")

fviz_pca_ind(pca_result1, habillage = groups, shape = 16, mean.point=F, 
             addEllipses = T, legend.title="Groups", ellipse.type="confidence", ellipse.level=0.99, 
             palette = c("CK" = "#f19a8e", "Y0M_2d" = "#e377c2", "Y0_25" = "#9cd9e8", 
                         "Y0_5_1d" = "#71cabc", "Y0_5_2d" = "#92a0bd", 
                         "Y0_5_4d" = "#f9c7b8", "Y0_75" = "#bac1d5"))+ #Cell配色哦 
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加个边框


ggbiplot(pca_result,obs.scale = 1,var.scale = 1, groups = groups,
         ellipse = T,var.axes = F) + scale_color_brewer(palette = "Set1") + 
  theme(legend.direction = 'horizontal',legend.position = 'top') + 
  geom_vline(xintercept = c(0), linetype = 6,color = "black") + 
  geom_hline(yintercept = c(0),linetype = 6,color = "black") + 
  theme_bw() + theme(panel.grid = element_line(colour = NA))




dds <- read.csv("transcript.count.matrix_all.csv",row.names = 1)
dds <- dds+1
dds <- as.matrix(dds)
colnames(dds)
coldata <- data.frame(row.names=colnames(dds), groups)
# 没有提前将count矩阵的数值转换成整数的话加上round()函数
dds <- DESeqDataSetFromMatrix(countData = round(dds), colData = coldata, 
                              design = ~ groups)
rld <- rlog(dds)
rld_period <- plotPCA(rld, intgroup=c("groups"),returnData =
                        F)
rld_period <- data.frame(rld_period)
library(ggrepel)
library(ggalt)
library(ggforce)

p1<- ggplot()+
  labs(x = "PC1 48%",y = "PC2 21%")+
  labs(title = "nutrtion")+
  geom_hline(yintercept =0,linetype =2,size =2)+
  geom_vline(xintercept = 0,linetype  = 2,size = 2)+
  theme(panel.grid = element_blank())+
  geom_point(data = rld_period,aes(x = PC1, y= PC2,col = groups),size = 5)+
  scale_shape_manual(values = c(15,19,17,18,25))+
  scale_fill_npg()+
  labs(color = "Time",shape = "Type")+
  stat_ellipse(data = rld_period,aes(x = PC1,y = PC2,group = period),level = 0.95,
               linetype = 2,lwd = 0.1,color =1)
#添加圈置信圈
mytheme<-theme_bw() + theme(axis.title = element_text(size = 20), 
                            axis.text = element_text(size = 20), 
                            panel.grid.major = element_line(color ="white"), 
                            panel.grid.minor = element_line(colour = "white"), 
                            axis.text.x = element_text(size = 15,angle = 0,
                                                       vjust = 0, hjust = 0.5,color = "black"),
                            axis.text.y = element_text(size = 15,color ="black"),
                            legend.text = element_text(size = 15),
                            #legend.key.size = unit(10,"pt"),
                            legend.title = element_blank(),
                            legend.direction = "vertical",
                            #legend.position = c(.85,.75),
                            panel.border = element_rect(fill = NA, color = "black", 
                                                        size = 3, linetype ="solid"),
                            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))########图像边距




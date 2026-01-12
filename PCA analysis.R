
library(ggplot2)
BiocManager::install("ggalt")
library(factoextra)
library(FactoMineR)
library(ggbiplot)

data_path <- "RNA_seq_1week/salt-stress_R/PCA/transcript.tpm.matrix_all.csv"
data <- read.csv(data_path, header = TRUE, row.names = 1)

data_t <- t(data)

data_t[is.na(data_t)] <- 0  
data_t[!is.finite(data_t)] <- 0  
data_t <- data_t[rowSums(data_t) > 0, ]
data_t <- data_t[, apply(data_t, 2, var) > 0]
data_log <- log2(data_t + 1)

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
    groups[i] <- name  
  }
}
groups <- factor(groups, levels = c("CK", "Y0M_2d", "Y0_25",  "Y0_5_1d", "Y0_5_2d", "Y0_5_4d", "Y0_75"))

pca_result1 <- prcomp(data_log, center = TRUE, scale. = TRUE)

fviz_pca_ind(pca_result1, habillage = groups, shape = 16, mean.point=F, 
             addEllipses = T, legend.title="Groups", ellipse.type="confidence", ellipse.level=0.99, 
             palette = c("CK" = "#f19a8e", "Y0M_2d" = "#e377c2", "Y0_25" = "#9cd9e8", 
                         "Y0_5_1d" = "#71cabc", "Y0_5_2d" = "#92a0bd", 
                         "Y0_5_4d" = "#f9c7b8", "Y0_75" = "#bac1d5"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))





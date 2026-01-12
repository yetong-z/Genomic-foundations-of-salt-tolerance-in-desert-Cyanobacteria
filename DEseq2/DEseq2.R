BiocManager::install('rstatix')
library('DESeq2')
library(tidyverse)
library(rstatix)


result_dir <- "Result"

FCcut <- 2
FDRcut <- 0.05

# get counts matrix
count_files <- list.files(pattern = "^count_matrix.*\\.csv$")

print(count_files)

if (length(count_files) == 0) {
  stop("no files")
}

for (count_file in count_files) {
  # get file name
  file_base <- tools::file_path_sans_ext(count_file)
  
  # contract comparison information
  comparison_info <- gsub("count_matrix_", "", file_base)
  
  cat("\nfile:", count_file, "\n")
  cat("comparison information:", comparison_info, "\n")
  
  dds <- read.csv(count_file, row.names = 1)
  
  # get conditions
  groups <- unlist(strsplit(comparison_info, "_"))
  
  if (length(groups) >= 2) {
    sample_count <- ncol(dds) / length(groups)
    
    if (sample_count %% 1 != 0) {
      warning("The sample size cannot be evenly distributed among each group")
      # Default packet deny
      sample_count <- floor(ncol(dds) / 2)
      count_condition <- factor(c(rep(groups[1], sample_count), 
                                rep(paste(groups[-1], collapse="_"), ncol(dds) - sample_count)))
    } else {
      # create condition factors
      count_condition <- factor(unlist(lapply(1:length(groups), function(i) {
        rep(groups[i], sample_count)
      })))
    }
  } else {
    warning("The comparison group information cannot be parsed")
    count_condition <- factor(c(rep("Group1", floor(ncol(dds)/2)), 
                              rep("Group2", ncol(dds) - floor(ncol(dds)/2))))
  }
  
  cat("groups:\n")
  print(count_condition)
  
  coldata <- data.frame(row.names=colnames(dds), count_condition)
  
  # create database
  dds_obj <- DESeqDataSetFromMatrix(countData = round(dds), colData = coldata, 
                                design = ~ count_condition)
  
  dds_obj <- dds_obj[rowSums(counts(dds_obj)) > 1, ]
  
  # calculate sizefactors
  dds_obj.sizefactor <- estimateSizeFactors(dds_obj) 
  
  cat("Size factors:\n")
  print(sizeFactors(dds_obj.sizefactor))
  
  deq <- DESeq(dds_obj)
  res <- results(deq)
  
  resdata <- merge(as.data.frame(res), as.data.frame(counts(deq, normalized=TRUE)),
                  by="row.names", sort=FALSE)
  
  # save results
  output_file <- file.path(result_dir, paste0("deseq2_", comparison_info, ".csv"))
  write.csv(x = resdata, file = output_file)
  cat("save:", output_file, "\n")
  
  # contract DEGs
  out <- cbind(res$log2FoldChange, res$pvalue, res$padj)
  rownames(out) <- row.names(res)
  colnames(out) <- c('log2(FoldChange)', 'Pvalue', 'FDR')
  
  diff <- out[(!is.na(res$padj) & res$padj < FDRcut) & 
              abs(res$log2FoldChange) > abs(log2(FCcut)), ]
  
  all_result_file <- file.path(result_dir, 
                              paste0("DESeq2_AllResult_", comparison_info, ".csv"))
  write.table(out, all_result_file, sep = ",", quote=F, 
              row.names=T, col.names=T)
  cat("save all results:", all_result_file, "\n")
  
  deg_file <- file.path(result_dir, 
                       paste0("DESeq2_DEG_", comparison_info, ".csv"))
  write.table(diff, deg_file, sep = ",", 
              row.names=TRUE, quote=F)
  cat("save DEGs:", deg_file, "\n")
  
  deg_list_file <- file.path(result_dir, 
                            paste0("DESeq2_DEGlist_", comparison_info, ".csv"))
  write.table(rownames(diff), deg_list_file, sep = ",", 
              quote=F, row.names=F, col.names=F)
  cat("save DEGs expression list:", deg_list_file, "\n")
  
}

cat("\nDone！Results save in", result_dir, "\n")



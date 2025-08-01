library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

########################################

res_list <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_results_for_mashr.rds")

########################################

corr_df <- data.frame()

for (test in names(res_list)){
  df <- res_list[[test]]
  
  corr_df <- rbind(corr_df, data.frame(comparison = test, 
                                       spearman_beta = as.numeric(cor(df[2], df[6], method = "spearman")), 
                                       spearman_zscore = as.numeric(cor(df[11], df[12], method = "spearman"))))
}

saveRDS(corr_df, file = "./analysis/oligodendrocytes/differential_expression/DESeq_nebula_correlation.rds")

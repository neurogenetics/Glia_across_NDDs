library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

# read in merged DE results files

nebula_results <- list.files("./analysis/oligodendrocytes/differential_expression/nebula_results/", pattern = "\\.rds$", full.names = T)
nebula_results <- lapply(nebula_results, readRDS)
nebula_names <- list.files("./analysis/oligodendrocytes/differential_expression/nebula_results/", pattern = "\\.rds$", full.names = F)
nebula_names <- gsub("_nebula_results.rds", "", nebula_names)
names(nebula_results) <- nebula_names


DESeq_results <- list.files("./analysis/oligodendrocytes/differential_expression/DESeq_results/", pattern = "\\.rds$", full.names = T)
DESeq_results <- lapply(DESeq_results, readRDS)
DESeq_names <- list.files("./analysis/oligodendrocytes/differential_expression/DESeq_results/", pattern = "\\.rds$", full.names = F)
DESeq_names <- gsub("_DESeq_results.rds", "", DESeq_names)
names(DESeq_results) <- DESeq_names

##################################################

# organize / unify DESeq / nebula results

comparisons <- names(nebula_results)
# comparison = comparisons[1]

merged_results_list <- list()

for (comparison in comparisons) {
  # filter nebula results for genes w/ convergence = 1
  nebula_summary <- nebula_results[[comparison]]$summary
  nebula_summary$convergence <- nebula_results[[comparison]]$convergence
  nebula_summary <- nebula_summary[nebula_summary$convergence == 1, ]
  disease <- sub(".*_(.*)", "\\1", comparison)
  nebula_summary <- nebula_summary[, c("gene", 
                                       paste0("logFC_group", disease),
                                       paste0("se_group", disease),
                                       paste0("p_group", disease))]
  
  # filter for nebula SE < 10
  nebula_summary <- nebula_summary[nebula_summary[3] < 10, ]
  
  
  DESeq_res <- DESeq_results[[comparison]]
  DESeq_res <- na.omit(DESeq_res)
  
  # filter DESeq and nebula results for the same genes
  DESeq_res <- DESeq_res[DESeq_res$gene %in% nebula_summary$gene, ]
  nebula_summary <- nebula_summary[nebula_summary$gene %in% DESeq_res$gene, ]
  
  merged_res <- nebula_summary %>%
    left_join(DESeq_res, by = "gene")
  
  merged_res$nebula_zscore <- scale(merged_res[,2])
  merged_res$DESeq_zscore <- scale(merged_res[,6])
  
  merged_results_list[[comparison]] <- merged_res
}

saveRDS(merged_results_list, file = "./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_results_for_mashr.rds")

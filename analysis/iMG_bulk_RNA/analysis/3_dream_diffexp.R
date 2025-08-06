# done in interactive job w/ params for multithreading
# sinteractive --mem=500g --gres=lscratch:5 --cpus-per-task=24 --tunnel --time=2160

library(variancePartition)
library(tidyverse)
library(edgeR)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

cts <- readRDS("./analysis/cellculture/iMG_bulk_RNA/cleaned_counts/all_batches_cts_merged.rds")
meta <- readRDS("./analysis/cellculture/iMG_bulk_RNA/metadata/sample_meta_merged.rds")

##################################################
##################################################
##################################################

# run dream on different batches in a loop

# parallel processing 
param <- SnowParam(12, "SOCK", progressbar = TRUE)

batch_list <- c("XR10585", "DA10634")
# batch_list <- c("DA10634")

##################################################

# run WITHOUT covariates

for (batch in batch_list){
  meta_filtered <- meta[grepl(batch, rownames(meta)), ]
  cts_filtered <- cts[,colnames(cts) %in% rownames(meta_filtered)]
  
  
  # filter for expressed genes
  # filtering so genes are expressed @ cpm > 1 in at least 3 samples of any treatment group
  dge <- DGEList(cts_filtered)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm <- as.data.frame(edgeR::cpm(dge))
  
  meta_filtered <- meta_filtered[match(colnames(cpm), rownames(meta_filtered)), ]
  
  genes_keep <- vector("logical", nrow(cpm))
  genes_per_group <- data.frame(treatment = character(), n_genes = numeric(), stringsAsFactors = FALSE)
  
  for (t in unique(meta_filtered$treatment)) {
    samples_in_group <- rownames(meta_filtered)[meta_filtered$treatment == t]
    
    cpm_sub <- cpm[, samples_in_group, drop = FALSE]
    
    gene_pass <- rowSums(cpm_sub > 1) >= 3
    
    genes_keep <- genes_keep | gene_pass
    
    genes_per_group <- rbind(genes_per_group,
                             data.frame(treatment = t, n_genes = sum(gene_pass), stringsAsFactors = FALSE))
  }
  
  gene_pass <- rownames(cpm)[gene_pass]
  
  cts_filtered <- cts_filtered[rownames(cts_filtered) %in% gene_pass, ]
  
  
  # loop to run dream
  treatments <- unique(meta_filtered$treatment_dose)
  treatments <- treatments[!treatments %in% "CTL_0"]
  
  res_list <- list()
  res_obj_list <- list()
  
  for (treatment in treatments){
    print(paste0("Processing control vs. ", treatment))
    
    # filter for the samples in this comparison
    meta_test <- meta_filtered[meta_filtered$treatment_dose %in% c("CTL_0", treatment), ]
    cts_test <- cts_filtered[, colnames(cts_filtered) %in% rownames(meta_test)]
    meta_test$PPMI_line <- as.factor(meta_test$PPMI_line)
    meta_test$treatment_dose <- factor(meta_test$treatment_dose, levels = c("CTL_0", treatment))
    
    # dream workflow
    dge <- DGEList(cts_test)
    dge <- calcNormFactors(dge, method = "TMM")
    
    form <- ~ treatment_dose + (1 | PPMI_line)
    
    vobjDream <- voomWithDreamWeights(dge, form, meta_test, BPPARAM = param)
    fitmm <- dream(vobjDream, form, meta_test, BPPARAM = param)
    fitmm <- eBayes(fitmm)
    
    res <- topTable(fitmm, coef = paste0("treatment_dose", treatment), number = Inf)
    
    res_list[[treatment]] <- res
    res_obj_list[[treatment]] <- fitmm
  }
  saveRDS(res_list, paste0("./analysis/cellculture/iMG_bulk_RNA/tmp_results/", batch, "_dream_res_df_list.rds"))
  saveRDS(res_obj_list, paste0("./analysis/cellculture/iMG_bulk_RNA/tmp_results/", batch, "_dream_res_obj_list.rds"))
}

##################################################

# run WITH covariates (% intronic, % ribo, % mito)

for (batch in batch_list){
  meta_filtered <- meta[grepl(batch, rownames(meta)), ]
  cts_filtered <- cts[,colnames(cts) %in% rownames(meta_filtered)]
  
  
  # filter for expressed genes
  # filtering so genes are expressed @ cpm > 1 in at least 3 samples of any treatment group
  dge <- DGEList(cts_filtered)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm <- as.data.frame(edgeR::cpm(dge))
  
  meta_filtered <- meta_filtered[match(colnames(cpm), rownames(meta_filtered)), ]
  
  genes_keep <- vector("logical", nrow(cpm))
  genes_per_group <- data.frame(treatment = character(), n_genes = numeric(), stringsAsFactors = FALSE)
  
  for (t in unique(meta_filtered$treatment)) {
    samples_in_group <- rownames(meta_filtered)[meta_filtered$treatment == t]
    
    cpm_sub <- cpm[, samples_in_group, drop = FALSE]
    
    gene_pass <- rowSums(cpm_sub > 1) >= 3
    
    genes_keep <- genes_keep | gene_pass
    
    genes_per_group <- rbind(genes_per_group,
                             data.frame(treatment = t, n_genes = sum(gene_pass), stringsAsFactors = FALSE))
  }
  
  gene_pass <- rownames(cpm)[gene_pass]
  
  cts_filtered <- cts_filtered[rownames(cts_filtered) %in% gene_pass, ]
  
  
  # loop to run dream
  treatments <- unique(meta_filtered$treatment_dose)
  treatments <- treatments[!treatments %in% "CTL_0"]
  
  res_list <- list()
  res_obj_list <- list()
  
  for (treatment in treatments){
    print(paste0("Processing control vs. ", treatment))
    
    # filter for the samples in this comparison
    meta_test <- meta_filtered[meta_filtered$treatment_dose %in% c("CTL_0", treatment), ]
    cts_test <- cts_filtered[, colnames(cts_filtered) %in% rownames(meta_test)]
    meta_test$PPMI_line <- as.factor(meta_test$PPMI_line)
    meta_test$treatment_dose <- factor(meta_test$treatment_dose, levels = c("CTL_0", treatment))
    meta_test$sex <- as.factor(meta_test$sex)
    # meta_test$nanodrop_scaled <- scale(meta_test$RNA_nanodrop)
    meta_test$intronic_scaled <- scale(meta_test$pct.intronic)
    meta_test$mito_scaled <- scale(meta_test$pct.mito)
    meta_test$ribo_scaled <- scale(meta_test$pct.ribo)
    
    # dream workflow
    dge <- DGEList(cts_test)
    dge <- calcNormFactors(dge, method = "TMM")
    
    # including [RNA]nanodrop (either scaled or unscaled) causes immediate error
    # including sex as either random or fixed effect causes error after iterations (maybe bc co-linear w/ line?)
    form <- ~ treatment_dose + intronic_scaled + mito_scaled + ribo_scaled + (1 | PPMI_line)
    
    vobjDream <- voomWithDreamWeights(dge, form, meta_test, BPPARAM = param)
    fitmm <- dream(vobjDream, form, meta_test, BPPARAM = param)
    fitmm <- eBayes(fitmm)
    
    res <- topTable(fitmm, coef = paste0("treatment_dose", treatment), number = Inf)
    
    res_list[[treatment]] <- res
    res_obj_list[[treatment]] <- fitmm
  }
  saveRDS(res_list, paste0("./analysis/cellculture/iMG_bulk_RNA/tmp_results/", batch, "_dream_res_df_list_w_covars.rds"))
  saveRDS(res_obj_list, paste0("./analysis/cellculture/iMG_bulk_RNA/tmp_results/", batch, "_dream_res_obj_list_w_covars.rds"))
}

library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################
##################################################
##################################################

# getting gene scores for each GEP and top 200 genes for each GEP for each study

##################################################

# AMP-PD

gene_scores_AMPPD <- read.delim("./analysis/microglia/cluster_characterization/cNMF_by_dataset/AMPPD/microglia/microglia.gene_spectra_score.k_17.dt_0_12.txt", header = T)
gene_scores_AMPPD <- as.data.frame(t(gene_scores_AMPPD))
gene_scores_AMPPD <- gene_scores_AMPPD[!rownames(gene_scores_AMPPD) == "X", ]
colnames(gene_scores_AMPPD) <- paste0("GEP", 1:17)
gene_scores_AMPPD <- na.omit(gene_scores_AMPPD)

# get top 200 genes from each GEP
GEP_list_AMPPD <- lapply(1:ncol(gene_scores_AMPPD), function(i) {
  top_indices <- order(gene_scores_AMPPD[, i], decreasing = TRUE)[1:200]
  data.frame(gene = rownames(gene_scores_AMPPD)[top_indices], enrichment = gene_scores_AMPPD[top_indices, i])
})
names(GEP_list_AMPPD) <- colnames(gene_scores_AMPPD)

saveRDS(GEP_list_AMPPD, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_top200_list.rds")
saveRDS(gene_scores_AMPPD, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_gene_scores.rds")

##################################################

# Gerrits

gene_scores_Gerrits <- read.delim("./analysis/microglia/cluster_characterization/cNMF_by_dataset/Gerrits/microglia/microglia.gene_spectra_score.k_16.dt_0_16.txt", header = T)
gene_scores_Gerrits <- as.data.frame(t(gene_scores_Gerrits))
gene_scores_Gerrits <- gene_scores_Gerrits[!rownames(gene_scores_Gerrits) == "X", ]
colnames(gene_scores_Gerrits) <- paste0("GEP", 1:16)
gene_scores_Gerrits <- na.omit(gene_scores_Gerrits)

# get top 200 genes from each GEP
GEP_list_Gerrits <- lapply(1:ncol(gene_scores_Gerrits), function(i) {
  top_indices <- order(gene_scores_Gerrits[, i], decreasing = TRUE)[1:200]
  data.frame(gene = rownames(gene_scores_Gerrits)[top_indices], enrichment = gene_scores_Gerrits[top_indices, i])
})
names(GEP_list_Gerrits) <- colnames(gene_scores_Gerrits)

saveRDS(GEP_list_Gerrits, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_top200_list.rds")
saveRDS(gene_scores_Gerrits, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_gene_scores.rds")

##################################################

# Mathys

gene_scores_Mathys <- read.delim("./analysis/microglia/cluster_characterization/cNMF_by_dataset/Mathys/microglia/microglia.gene_spectra_score.k_12.dt_0_12.txt", header = T)
gene_scores_Mathys <- as.data.frame(t(gene_scores_Mathys))
gene_scores_Mathys <- gene_scores_Mathys[!rownames(gene_scores_Mathys) == "X", ]
colnames(gene_scores_Mathys) <- paste0("GEP", 1:12)
gene_scores_Mathys <- na.omit(gene_scores_Mathys)

# get top 200 genes from each GEP
GEP_list_Mathys <- lapply(1:ncol(gene_scores_Mathys), function(i) {
  top_indices <- order(gene_scores_Mathys[, i], decreasing = TRUE)[1:200]
  data.frame(gene = rownames(gene_scores_Mathys)[top_indices], enrichment = gene_scores_Mathys[top_indices, i])
})
names(GEP_list_Mathys) <- colnames(gene_scores_Mathys)

saveRDS(GEP_list_Mathys, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_top200_list.rds")
saveRDS(gene_scores_Mathys, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_gene_scores.rds")

##################################################

# Pineda

gene_scores_Pineda <- read.delim("./analysis/microglia/cluster_characterization/cNMF_by_dataset/Pineda/microglia/microglia.gene_spectra_score.k_12.dt_0_1.txt", header = T)
gene_scores_Pineda <- as.data.frame(t(gene_scores_Pineda))
gene_scores_Pineda <- gene_scores_Pineda[!rownames(gene_scores_Pineda) == "X", ]
colnames(gene_scores_Pineda) <- paste0("GEP", 1:12)
gene_scores_Pineda <- na.omit(gene_scores_Pineda)

# get top 200 genes from each GEP
GEP_list_Pineda <- lapply(1:ncol(gene_scores_Pineda), function(i) {
  top_indices <- order(gene_scores_Pineda[, i], decreasing = TRUE)[1:200]
  data.frame(gene = rownames(gene_scores_Pineda)[top_indices], enrichment = gene_scores_Pineda[top_indices, i])
})
names(GEP_list_Pineda) <- colnames(gene_scores_Pineda)

saveRDS(GEP_list_Pineda, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_top200_list.rds")
saveRDS(gene_scores_Pineda, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_gene_scores.rds")

##################################################
##################################################
##################################################

# getting average GEP usage across clusters

##################################################

# AMP-PD

usage_AMPPD <- read.delim("./analysis/microglia/cluster_characterization/cNMF_by_dataset/AMPPD/microglia/microglia.usages.k_17.dt_0_12.consensus.txt", header = T)
usage_AMPPD <- usage_AMPPD %>% 
  dplyr::rename(barcode = X)
colnames(usage_AMPPD)[2:18] <- paste0("GEP", 1:17)

cell_meta_AMPPD <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")
cell_meta_AMPPD <- cell_meta_AMPPD[rownames(cell_meta_AMPPD) %in% usage_AMPPD$barcode, ]
cell_meta_AMPPD <- cell_meta_AMPPD %>%
  rownames_to_column(var = "barcode")
cell_meta_AMPPD <- cell_meta_AMPPD[c("barcode", "cluster_anno")]

usage_AMPPD <- usage_AMPPD %>%
  left_join(cell_meta_AMPPD, by = "barcode")

usage_avg_AMPPD <- usage_AMPPD %>%
  group_by(cluster_anno) %>%
  summarise(across(starts_with("GEP"), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "cluster_anno")

saveRDS(usage_avg_AMPPD, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_avg_cluster_usage.rds")

##################################################

# Gerrits

usage_Gerrits <- read.delim("./analysis/microglia/cluster_characterization/cNMF_by_dataset/Gerrits/microglia/microglia.usages.k_16.dt_0_16.consensus.txt", header = T)
usage_Gerrits <- usage_Gerrits %>% 
  dplyr::rename(barcode = X)
colnames(usage_Gerrits)[2:17] <- paste0("GEP", 1:16)

cell_meta_Gerrits <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")
cell_meta_Gerrits <- cell_meta_Gerrits[rownames(cell_meta_Gerrits) %in% usage_Gerrits$barcode, ]
cell_meta_Gerrits <- cell_meta_Gerrits %>%
  rownames_to_column(var = "barcode")
cell_meta_Gerrits <- cell_meta_Gerrits[c("barcode", "cluster_anno")]

usage_Gerrits <- usage_Gerrits %>%
  left_join(cell_meta_Gerrits, by = "barcode")

usage_avg_Gerrits <- usage_Gerrits %>%
  group_by(cluster_anno) %>%
  summarise(across(starts_with("GEP"), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "cluster_anno")

saveRDS(usage_avg_Gerrits, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_avg_cluster_usage.rds")

##################################################

# Mathys

usage_Mathys <- read.delim("./analysis/microglia/cluster_characterization/cNMF_by_dataset/Mathys/microglia/microglia.usages.k_12.dt_0_12.consensus.txt", header = T)
usage_Mathys <- usage_Mathys %>% 
  dplyr::rename(barcode = X)
colnames(usage_Mathys)[2:13] <- paste0("GEP", 1:12)

cell_meta_Mathys <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")
cell_meta_Mathys <- cell_meta_Mathys[rownames(cell_meta_Mathys) %in% usage_Mathys$barcode, ]
cell_meta_Mathys <- cell_meta_Mathys %>%
  rownames_to_column(var = "barcode")
cell_meta_Mathys <- cell_meta_Mathys[c("barcode", "cluster_anno")]

usage_Mathys <- usage_Mathys %>%
  left_join(cell_meta_Mathys, by = "barcode")

usage_avg_Mathys <- usage_Mathys %>%
  group_by(cluster_anno) %>%
  summarise(across(starts_with("GEP"), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "cluster_anno")

saveRDS(usage_avg_Mathys, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_avg_cluster_usage.rds")

##################################################

# Pineda

usage_Pineda <- read.delim("./analysis/microglia/cluster_characterization/cNMF_by_dataset/Pineda/microglia/microglia.usages.k_12.dt_0_1.consensus.txt", header = T)
usage_Pineda <- usage_Pineda %>% 
  dplyr::rename(barcode = X)
colnames(usage_Pineda)[2:13] <- paste0("GEP", 1:12)

cell_meta_Pineda <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")
cell_meta_Pineda <- cell_meta_Pineda[rownames(cell_meta_Pineda) %in% usage_Pineda$barcode, ]
cell_meta_Pineda <- cell_meta_Pineda %>%
  rownames_to_column(var = "barcode")
cell_meta_Pineda <- cell_meta_Pineda[c("barcode", "cluster_anno")]

usage_Pineda <- usage_Pineda %>%
  left_join(cell_meta_Pineda, by = "barcode")

usage_avg_Pineda <- usage_Pineda %>%
  group_by(cluster_anno) %>%
  summarise(across(starts_with("GEP"), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "cluster_anno")

saveRDS(usage_avg_Pineda, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_avg_cluster_usage.rds")

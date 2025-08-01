library(Seurat)
library(tidyverse)
library(scCustomize)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs")

########################################

# function for processing
source("./code_organized/functions/seurat_lognorm_harmony_integration.R")

########################################
########################################
########################################

# AMP-PD

########################################

PD_DMV_oligos <- readRDS("./oligos_data/subset_by_dataset_region/PD_DMV_oligos.rds")
PD_GPi_oligos <- readRDS("./oligos_data/subset_by_dataset_region/PD_GPi_oligos.rds")
PD_M1_oligos <- readRDS("./oligos_data/subset_by_dataset_region/PD_M1_oligos.rds")
PD_PFC_oligos <- readRDS("./oligos_data/subset_by_dataset_region/PD_PFC_oligos.rds")
PD_V1_oligos <- readRDS("./oligos_data/subset_by_dataset_region/PD_V1_oligos.rds")

########################################

# need to process these separately bc there are too many cells

########################################

# DMV
PD_DMV_oligos <- lognorm_harmonyintegrate(PD_DMV_oligos, n_varfeats = 2000, vars_to_regress = "pct.mito", n_dims = NULL, 
                                          cluster_res = 0.1, harmony_batches = c("donor", "batch"))

PD_DMV_oligos_cts <- LayerData(PD_DMV_oligos, layer = "counts", assay = "RNA")
PD_DMV_oligos_meta <- PD_DMV_oligos@meta.data

saveRDS(PD_DMV_oligos_cts, file = "./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_DMV_oligos_cts.rds")
saveRDS(PD_DMV_oligos_meta, file = "./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_DMV_oligos_meta.rds")

########################################

# GPi
PD_GPi_oligos <- lognorm_harmonyintegrate(PD_GPi_oligos, n_varfeats = 2000, vars_to_regress = "pct.mito", n_dims = NULL, 
                                          cluster_res = 0.1, harmony_batches = c("donor", "batch"))

PD_GPi_oligos_cts <- LayerData(PD_GPi_oligos, layer = "counts", assay = "RNA")
PD_GPi_oligos_meta <- PD_GPi_oligos@meta.data

saveRDS(PD_GPi_oligos_cts, file = "./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_GPi_oligos_cts.rds")
saveRDS(PD_GPi_oligos_meta, file = "./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_GPi_oligos_meta.rds")

########################################

# M1
PD_M1_oligos <- lognorm_harmonyintegrate(PD_M1_oligos, n_varfeats = 2000, vars_to_regress = "pct.mito", n_dims = NULL,
                                         cluster_res = 0.1, harmony_batches = c("donor", "batch"))

PD_M1_oligos_cts <- LayerData(PD_M1_oligos, layer = "counts", assay = "RNA")
PD_M1_oligos_meta <- PD_M1_oligos@meta.data

saveRDS(PD_M1_oligos_cts, file = "./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_M1_oligos_cts.rds")
saveRDS(PD_M1_oligos_meta, file = "./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_M1_oligos_meta.rds")

########################################

# PFC
PD_PFC_oligos <- lognorm_harmonyintegrate(PD_PFC_oligos, n_varfeats = 2000, vars_to_regress = "pct.mito", n_dims = NULL,
                                         cluster_res = 0.1, harmony_batches = c("donor", "batch"))

PD_PFC_oligos_cts <- LayerData(PD_PFC_oligos, layer = "counts", assay = "RNA")
PD_PFC_oligos_meta <- PD_PFC_oligos@meta.data

saveRDS(PD_PFC_oligos_cts, file = "./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_PFC_oligos_cts.rds")
saveRDS(PD_PFC_oligos_meta, file = "./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_PFC_oligos_meta.rds")

########################################

# V1
PD_V1_oligos <- lognorm_harmonyintegrate(PD_V1_oligos, n_varfeats = 2000, vars_to_regress = "pct.mito", n_dims = NULL,
                                         cluster_res = 0.1, harmony_batches = c("donor", "batch"))

PD_V1_oligos_cts <- LayerData(PD_V1_oligos, layer = "counts", assay = "RNA")
PD_V1_oligos_meta <- PD_V1_oligos@meta.data

saveRDS(PD_V1_oligos_cts, file = "./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_V1_oligos_cts.rds")
saveRDS(PD_V1_oligos_meta, file = "./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_V1_oligos_meta.rds")

########################################
########################################
########################################

# Mathys

########################################

AD_AnG_oligos <- readRDS("./oligos_data/subset_by_dataset_region/AD_AnG_oligos.rds")
AD_EC_oligos <- readRDS("./oligos_data/subset_by_dataset_region/AD_EC_oligos.rds")
AD_HIP_oligos <- readRDS("./oligos_data/subset_by_dataset_region/AD_HIP_oligos.rds")
AD_MTG_oligos <- readRDS("./oligos_data/subset_by_dataset_region/AD_MTG_oligos.rds")
AD_PFC_oligos <- readRDS("./oligos_data/subset_by_dataset_region/AD_PFC_oligos.rds")
AD_TH_oligos <- readRDS("./oligos_data/subset_by_dataset_region/AD_TH_oligos.rds")

########################################

AD_oligos_list <- list(AD_AnG_oligos, AD_EC_oligos, AD_HIP_oligos, AD_MTG_oligos, AD_PFC_oligos, AD_TH_oligos)
AD_oligos_merged <- Merge_Seurat_List(AD_oligos_list)

AD_oligos_merged <- lognorm_harmonyintegrate(AD_oligos_merged, n_varfeats = 2000, vars_to_regress = "pct.mito", n_dims = NULL,
                                             cluster_res = 0.1, harmony_batches = c("donor"))

AD_oligos_merged_cts <- LayerData(AD_oligos_merged, layer = "counts", assay = "RNA")
AD_oligos_merged_meta <- AD_oligos_merged@meta.data

saveRDS(AD_oligos_merged_cts, file = "./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AD_oligos_merged_cts.rds")
saveRDS(AD_oligos_merged_meta, file = "./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AD_oligos_merged_meta.rds")

########################################
########################################
########################################

# Pineda

########################################

ALS_M1_oligos <- readRDS("./oligos_data/subset_by_dataset_region/ALS_M1_oligos.rds")
ALS_PFC_oligos <- readRDS("./oligos_data/subset_by_dataset_region/ALS_PFC_oligos.rds")

########################################

ALS_oligos_list <- list(ALS_M1_oligos, ALS_PFC_oligos)
ALS_oligos_merged <- Merge_Seurat_List(ALS_oligos_list)

ALS_oligos_merged <- lognorm_harmonyintegrate(ALS_oligos_merged, n_varfeats = 2000, vars_to_regress = "pct.mito", n_dims = NULL,
                                             cluster_res = 0.1, harmony_batches = c("donor"))

ALS_oligos_merged_cts <- LayerData(ALS_oligos_merged, layer = "counts", assay = "RNA")
ALS_oligos_merged_meta <- ALS_oligos_merged@meta.data

saveRDS(ALS_oligos_merged_cts, file = "./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/ALS_oligos_merged_cts.rds")
saveRDS(ALS_oligos_merged_meta, file = "./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/ALS_oligos_merged_meta.rds")

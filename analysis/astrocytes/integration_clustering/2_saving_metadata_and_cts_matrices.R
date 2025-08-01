library(tidyverse)
library(Seurat)
library(scCustomize)
library(BPCells)

setwd("/data/ADRD/glia_across_NDDs")

########################################

astro <- readRDS("./combined_data/final_objects/astro_sketched_clustered_projected_NO_ctsfilter_res_0.15.rds")

astro_celllevel_meta <- astro@meta.data %>%
  dplyr::select(-orig.ident) %>%
  dplyr::select(-dataset_region)

# edit region names in metadata
astro_celllevel_meta <- astro_celllevel_meta %>%
  mutate(region = case_when(
    region == "PFC" ~ "FC/PFC",
    region == "FC" ~ "FC/PFC",
    region == "TC" ~ "TC/MTG",
    region == "MTG" ~ "TC/MTG",
    region == "OC" ~ "OC/V1", 
    region == "V1" ~ "OC/V1", 
    region == "HC" ~ "HIP",
    TRUE ~ region))

astro_celllevel_meta <- astro_celllevel_meta[!astro_celllevel_meta$celltype.full == "7", ]

saveRDS(astro_celllevel_meta, file = "./analysis/astrocytes/all_astrocytes_celllevel_metadata.rds")

########################################

# extract and save single-cell counts matrices from the different datasets

PD_astrocytes <- readRDS("./combined_data/processed_celltypes_by_dataset/PD_astrocytes_allregions_filtered.rds")
AD_astrocytes <- readRDS("./combined_data/processed_celltypes_by_dataset/AD_astrocytes_allregions_filtered.rds")
ALS_astrocytes <- readRDS("./combined_data/processed_celltypes_by_dataset/ALS_astrocytes_allregions_filtered.rds")
FTD_astrocytes <- readRDS("./combined_data/processed_celltypes_by_dataset/FTD_astrocytes_allregions_filtered.rds")

PD_counts <- GetAssayData(PD_astrocytes, assay = "RNA", layer = "counts")
AD_counts <- GetAssayData(AD_astrocytes, assay = "RNA", layer = "counts")
ALS_counts <- GetAssayData(ALS_astrocytes, assay = "RNA", layer = "counts")
FTD_counts <- GetAssayData(FTD_astrocytes, assay = "RNA", layer = "counts")

# filter out cells from cluster 7 (see above)
PD_counts_filtered <- PD_counts[, colnames(PD_counts) %in% rownames(astro_celllevel_meta)]
AD_counts_filtered <- AD_counts[, colnames(AD_counts) %in% rownames(astro_celllevel_meta)]
ALS_counts_filtered <- ALS_counts[, colnames(ALS_counts) %in% rownames(astro_celllevel_meta)]
FTD_counts_filtered <- FTD_counts[, colnames(FTD_counts) %in% rownames(astro_celllevel_meta)]

saveRDS(PD_counts, file = "./analysis/astrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_astrocytes_sc_cts.rds")
saveRDS(AD_counts, file = "./analysis/astrocytes/differential_expression/single_cell_cts_matrices/Mathys_astrocytes_sc_cts.rds")
saveRDS(ALS_counts, file = "./analysis/astrocytes/differential_expression/single_cell_cts_matrices/Pineda_astrocytes_sc_cts.rds")
saveRDS(FTD_counts, file = "./analysis/astrocytes/differential_expression/single_cell_cts_matrices/Gerrits_astrocytes_sc_cts.rds")

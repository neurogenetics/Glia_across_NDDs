library(Seurat)
library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs/combined_data")

########################################

# PD

PD_donor_meta <- readRDS("/data/ADRD/glia_across_NDDs/metadata/cleaned/PD_donor_meta_cleaned.rds")
PD_donor_meta <- PD_donor_meta[c("donor", "group")]

PD_micro <- readRDS("./processed_celltypes_by_dataset/PD_microglia_allregions_filtered.rds")
PD_astro <- readRDS("./processed_celltypes_by_dataset/PD_astrocytes_allregions_filtered.rds")

PD_micro_meta <- PD_micro@meta.data
PD_astro_meta <- PD_astro@meta.data

PD_micro_meta_merged <- PD_micro_meta %>%
  left_join(PD_donor_meta, by = "donor")
rownames(PD_micro_meta_merged) <- PD_micro_meta_merged$barcode_batch
rownames(PD_micro_meta_merged) <- paste0("__", rownames(PD_micro_meta_merged))
identical(rownames(PD_micro_meta_merged), rownames(PD_micro_meta))
PD_micro@meta.data <- PD_micro_meta_merged

PD_astro_meta_merged <- PD_astro_meta %>%
  left_join(PD_donor_meta, by = "donor")
rownames(PD_astro_meta_merged) <- PD_astro_meta_merged$barcode_batch
rownames(PD_astro_meta_merged) <- paste0("__", rownames(PD_astro_meta_merged))
identical(rownames(PD_astro_meta_merged), rownames(PD_astro_meta))
PD_astro@meta.data <- PD_astro_meta_merged

PD_micro_HC_PFC <- subset(PD_micro, subset = region == "PFC" & group == "HC")
PD_astro_HC_PFC <- subset(PD_astro, subset = region == "PFC" & group == "HC")

saveRDS(PD_micro_HC_PFC, "./HC_PFC_objects/PD_micro_HC_PFC.rds")
saveRDS(PD_astro_HC_PFC, "./HC_PFC_objects/PD_astro_HC_PFC.rds")

########################################

# ALS

ALS_donor_meta <- readRDS("/data/ADRD/glia_across_NDDs/metadata/cleaned/ALS_donor_meta_cleaned.rds")
ALS_donor_meta <- ALS_donor_meta[c("donor", "group")]
ALS_donor_meta$donor <- as.character(ALS_donor_meta$donor)

ALS_micro <- readRDS("./processed_celltypes_by_dataset/ALS_microglia_allregions_filtered.rds")
ALS_astro <- readRDS("./processed_celltypes_by_dataset/ALS_astrocytes_allregions_filtered.rds")

ALS_micro_meta <- ALS_micro@meta.data
ALS_astro_meta <- ALS_astro@meta.data

ALS_micro_meta_merged <- ALS_micro_meta %>%
  left_join(ALS_donor_meta, by = "donor")
rownames(ALS_micro_meta_merged) <- ALS_micro_meta_merged$barcode
rownames(ALS_micro_meta_merged) <- paste0("_", rownames(ALS_micro_meta_merged))
identical(rownames(ALS_micro_meta_merged), rownames(ALS_micro_meta))
ALS_micro@meta.data <- ALS_micro_meta_merged

ALS_astro_meta_merged <- ALS_astro_meta %>%
  left_join(ALS_donor_meta, by = "donor")
rownames(ALS_astro_meta_merged) <- ALS_astro_meta_merged$barcode
rownames(ALS_astro_meta_merged) <- paste0("_", rownames(ALS_astro_meta_merged))
identical(rownames(ALS_astro_meta_merged), rownames(ALS_astro_meta))
ALS_astro@meta.data <- ALS_astro_meta_merged

ALS_micro_HC_PFC <- subset(ALS_micro, subset = region == "PFC" & group == "PN")
ALS_astro_HC_PFC <- subset(ALS_astro, subset = region == "PFC" & group == "PN")

saveRDS(ALS_micro_HC_PFC, "./HC_PFC_objects/ALS_micro_HC_PFC.rds")
saveRDS(ALS_astro_HC_PFC, "./HC_PFC_objects/ALS_astro_HC_PFC.rds")

########################################

# AD

AD_donor_meta <- readRDS("/data/ADRD/glia_across_NDDs/metadata/cleaned/AD_donor_meta_cleaned.rds")
AD_donor_meta <- AD_donor_meta[c("subject", "group")]
AD_donor_meta$subject <- gsub("ROSMAP-", "", AD_donor_meta$subject)
AD_donor_meta <- AD_donor_meta %>%
  dplyr::rename(donor = subject)

AD_micro <- readRDS("./processed_celltypes_by_dataset/AD_microglia_allregions_filtered.rds")
AD_astro <- readRDS("./processed_celltypes_by_dataset/AD_astrocytes_allregions_filtered.rds")

AD_micro_meta <- AD_micro@meta.data
AD_astro_meta <- AD_astro@meta.data

AD_micro_meta_merged <- AD_micro_meta %>%
  left_join(AD_donor_meta, by = "donor")
rownames(AD_micro_meta_merged) <- AD_micro_meta_merged$barcode
rownames(AD_micro_meta_merged) <- paste0("__", rownames(AD_micro_meta_merged))
identical(rownames(AD_micro_meta_merged), rownames(AD_micro_meta))
AD_micro@meta.data <- AD_micro_meta_merged

AD_astro_meta_merged <- AD_astro_meta %>%
  left_join(AD_donor_meta, by = "donor")
rownames(AD_astro_meta_merged) <- AD_astro_meta_merged$barcode
rownames(AD_astro_meta_merged) <- paste0("__", rownames(AD_astro_meta_merged))
identical(rownames(AD_astro_meta_merged), rownames(AD_astro_meta))
AD_astro@meta.data <- AD_astro_meta_merged

AD_micro_HC_PFC <- subset(AD_micro, subset = region == "PFC" & group == "HC")
AD_astro_HC_PFC <- subset(AD_astro, subset = region == "PFC" & group == "HC")

saveRDS(AD_micro_HC_PFC, "./HC_PFC_objects/AD_micro_HC_PFC.rds")
saveRDS(AD_astro_HC_PFC, "./HC_PFC_objects/AD_astro_HC_PFC.rds")

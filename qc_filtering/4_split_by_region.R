library(Seurat)
library(tidyverse)
library(scCustomize)

########################################

# PD

pd_dir <- "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/pd"
pd_file_list <- list.files(path = pd_dir, pattern = "\\.rds$", full.names = TRUE)


DMV_list <- list()
GPi_list <- list()
M1_list <- list()
PFC_list <- list()
V1_list <- list()

regions <- c("DMNX", "GPI", "PMC", "PFC", "PVC")

# for testing: pd_file_list <- c("/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/pd/PM-PD_Set9_C2_cellbender_seurat_filtered.rds",
# "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/pd/PM-PD_Set84_E2_cellbender_seurat_filtered.rds")


for (file in pd_file_list){
  cellbender_seurat_filtered <- readRDS(file)
  
  regions_in_obj <- unique(cellbender_seurat_filtered$region)
  
  for (brainregion in regions){
    if (brainregion %in% regions_in_obj) {
      fastq_brainregion <- subset(cellbender_seurat_filtered, subset = region == brainregion)
      
      if (brainregion == "DMNX") {
        DMV_list[[length(DMV_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "GPI") {
        GPi_list[[length(GPi_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "PMC") {
        M1_list[[length(M1_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "PFC") {
        PFC_list[[length(PFC_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "PVC") {
        V1_list[[length(V1_list) + 1]] <- fastq_brainregion
      }
    }
  }
}

# merge regions into their own objects

DMV_merged <- Merge_Seurat_List(DMV_list)
saveRDS(DMV_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/DMV_merged.rds")

GPi_merged <- Merge_Seurat_List(GPi_list)
saveRDS(GPi_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/GPi_merged.rds")

M1_merged <- Merge_Seurat_List(M1_list)
saveRDS(M1_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/M1_merged.rds")

# need to remove this one because they only have 1 cell
V1_list <- V1_list[-53]
V1_list <- V1_list[-93]
V1_merged <- Merge_Seurat_List(V1_list)
saveRDS(V1_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/V1_merged.rds")

PFC_merged <- Merge_Seurat_List(PFC_list)
saveRDS(PFC_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/PFC_merged.rds")

rm(list = ls())

########################################

# ALS

als_dir <- "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/als"
als_file_list <- list.files(path = als_dir, pattern = "\\.rds$", full.names = TRUE)


M1_list <- list()
PFC_list <- list()

regions <- c("M1", "PFC")

# for testing: als_file_list <- c("/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/als/SRR27882219_cellbender_seurat_filtered.rds",
# "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/als/SRR27882078_cellbender_seurat_filtered.rds")

for (file in als_file_list){
  sample_id <- sub(".*(SRR\\d{8}).*", "\\1", basename(file))
  
  cellbender_seurat_filtered <- readRDS(file)
  
  cellbender_seurat_filtered$srr_id <- sample_id
  
  regions_in_obj <- unique(cellbender_seurat_filtered$region)
  
  for (brainregion in regions){
    if (brainregion %in% regions_in_obj) {
      fastq_brainregion <- subset(cellbender_seurat_filtered, subset = region == brainregion)
      
      if (brainregion == "M1") {
        M1_list[[length(M1_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "PFC") {
        PFC_list[[length(PFC_list) + 1]] <- fastq_brainregion
    }
   }
  }
}

colnames(M1_list[[43]]) <- paste0(colnames(M1_list[[43]]), "_2")
colnames(M1_list[[52]]) <- paste0(colnames(M1_list[[52]]), "_2")
M1_merged <- Merge_Seurat_List(M1_list)
saveRDS(M1_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/als/M1_merged.rds")

PFC_merged <- Merge_Seurat_List(PFC_list)
saveRDS(PFC_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/als/PFC_merged.rds")

rm(list = ls())

########################################

# AD

ad_dir <- "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ad"
ad_file_list <- list.files(path = ad_dir, pattern = "\\.rds$", full.names = TRUE)


AnG_list <- list()
TH_list <- list()
EC_list <- list()
MTG_list <- list()
HIP_list <- list()
PFC_list <- list()

regions <- c("AnG", "PFC", "TH", "EC", "HC", "MTG")

# for testing: ad_file_list <- c("/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ad/D19-12442_cellbender_seurat_filtered.rds",
# "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ad/D19-8375_cellbender_seurat_filtered.rds",
# "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ad/D19-1916_cellbender_seurat_filtered.rds",
# "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ad/D19-12416_cellbender_seurat_filtered.rds",
# "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ad/D19-12333_cellbender_seurat_filtered.rds",
# "/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ad/D19-12362_cellbender_seurat_filtered.rds")

for (file in ad_file_list){
  cellbender_seurat_filtered <- readRDS(file)
  
  regions_in_obj <- unique(cellbender_seurat_filtered$region)
  
  for (brainregion in regions){
    if (brainregion %in% regions_in_obj) {
      fastq_brainregion <- subset(cellbender_seurat_filtered, subset = region == brainregion)
      
      if (brainregion == "AnG") {
        AnG_list[[length(AnG_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "PFC") {
        PFC_list[[length(PFC_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "TH") {
        TH_list[[length(TH_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "EC") {
        EC_list[[length(EC_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "HC") {
        HIP_list[[length(HIP_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "MTG") {
        MTG_list[[length(MTG_list) + 1]] <- fastq_brainregion
      }
    }
  }
}

AnG_merged <- Merge_Seurat_List(AnG_list)
saveRDS(AnG_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/AnG_merged.rds")

PFC_merged <- Merge_Seurat_List(PFC_list)
saveRDS(PFC_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/PFC_merged.rds")

TH_merged <- Merge_Seurat_List(TH_list)
saveRDS(TH_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/TH_merged.rds")

EC_merged <- Merge_Seurat_List(EC_list)
saveRDS(EC_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/EC_merged.rds")

HIP_merged <- Merge_Seurat_List(HIP_list)
saveRDS(HIP_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/HIP_merged.rds")

MTG_merged <- Merge_Seurat_List(MTG_list)
saveRDS(MTG_merged, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/MTG_merged.rds")

rm(list = ls())

########################################

# gerrits FTD

FTD_meta_filtered <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/FTD_metadata_FILTERED.rds")
FTD_fastq_list <- unique(sub(".*_", "", rownames(FTD_meta_filtered)))

FTD_seurat_list <- list()

# fastq <- "C4A"

for (fastq in FTD_fastq_list){
  cellbender_seurat_dir <- paste0("/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ftd/", fastq, "_cellbender_seurat_filtered.rds")
  cellbender_seurat <- readRDS(cellbender_seurat_dir)
  
  FTD_seurat_list[[fastq]] <- cellbender_seurat
}

FTD_seurat_merged <- Merge_Seurat_List(FTD_seurat_list)

FTD_basic_meta <- read.csv("/data/ADRD/glia_across_NDDs/metadata/gerrits_ftd_basic_meta.csv", header = T)
FTD_basic_meta <- FTD_basic_meta %>%
  rename(sample = Sample)
FTD_meta <- FTD_seurat_merged@meta.data
FTD_meta$barcodes <- rownames(FTD_meta)

FTD_meta_merged <- FTD_meta %>%
  left_join(FTD_basic_meta, by = "sample") %>%
  column_to_rownames(var = "barcodes")

identical(rownames(FTD_meta), rownames(FTD_meta_merged))

FTD_seurat_merged@meta.data <- FTD_meta_merged


FTD_FC <- subset(FTD_seurat_merged, subset = Region == "FC")
FTD_OC <- subset(FTD_seurat_merged, subset = Region == "OC")
FTD_TC <- subset(FTD_seurat_merged, subset = Region == "TC")

saveRDS(FTD_FC, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ftd/FC_merged.rds")
saveRDS(FTD_OC, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ftd/OC_merged.rds")
saveRDS(FTD_TC, file = "/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ftd/TC_merged.rds")

library(Seurat)
library(tidyverse)

PD_meta_filtered <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/PD_metadata_FILTERED.rds")
ALS_meta_filtered <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/ALSFTD_metadata_FILTERED.rds")
AD_meta_filtered <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/AD_metadata_FILTERED.rds")
FTD_meta_filtered <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/FTD_metadata_FILTERED.rds")

PD_cells <- rownames(PD_meta_filtered)
ALS_cells <- rownames(ALS_meta_filtered)
AD_cells <- rownames(AD_meta_filtered)
FTD_cells <- rownames(FTD_meta_filtered)

########################################

# PD

PD_meta_filtered$batch <- gsub("PM-PD_", "", PD_meta_filtered$batch)
PD_meta_filtered$fastq <- paste0("PM-PD_", PD_meta_filtered$batch)
PD_fastq_list_filtered <- unique(PD_meta_filtered$fastq)

# filter CellBender seurat objects to only have cells to keep

# for testing: fastq <- "PM-PD_Set1_1"

for (fastq in PD_fastq_list_filtered){
  cellbender_seurat_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/", fastq, "/", fastq, "_cellbender_seurat.rds")
  cellbender_seurat <- readRDS(cellbender_seurat_dir)
  
  cells_in_sample <- Cells(cellbender_seurat)
  cells_to_keep_in_sample <- intersect(cells_in_sample, PD_cells)
  
  cellbender_seurat_filtered <- subset(cellbender_seurat, cells = cells_to_keep_in_sample)
  
  filtered_save_dir <- paste0("/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/pd/", fastq, "_cellbender_seurat_filtered.rds")
  saveRDS(cellbender_seurat_filtered, file = filtered_save_dir)
}

########################################

# ALS

ALS_fastq_list_filtered <- unique(ALS_meta_filtered$fastq)

# for testing: fastq <- "SRR27882127"

for (fastq in ALS_fastq_list_filtered){
  cellbender_seurat_path <- paste0("/data/ADRD/ALSFTD_multiregion/fastq_processing/final_outs/", fastq)
  cellbender_seurat_dir <- list.files(path = cellbender_seurat_path, pattern = "_cellbender_seurat.rds", full.names = TRUE)
  cellbender_seurat <- readRDS(cellbender_seurat_dir)
  
  cells_in_sample <- Cells(cellbender_seurat)
  cells_to_keep_in_sample <- intersect(cells_in_sample, ALS_cells)
  
  cellbender_seurat_filtered <- subset(cellbender_seurat, cells = cells_to_keep_in_sample)
  
  filtered_save_dir <- paste0("/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/als/", fastq, "_cellbender_seurat_filtered.rds")
  saveRDS(cellbender_seurat_filtered, file = filtered_save_dir)
}

########################################

# AD

AD_meta_filtered$sample_region <- paste0(AD_meta_filtered$sample, "_", AD_meta_filtered$region)
AD_sample_list_filtered <- unique(AD_meta_filtered$sample_region)

# checking to make sure there is representation from all the samples in the filtered data so for loop doesn't fail
# AD_basic_meta <- read.csv(file = "/data/ADRD/AD_multiregion/fastq_processing/AD_basic_metadata.csv", header = T)
# filtered list has all 282 samples

AD_fastq_list <- read.csv(file = "/data/ADRD/AD_multiregion/fastq_processing/round_1_fastqs.txt", sep = "\t", header = F)
AD_fastq_list <- AD_fastq_list$V1

# for testing: fastq <- "D17-9566"

for (fastq in AD_fastq_list){
  cellbender_seurat_path <- paste0("/data/ADRD/AD_multiregion/fastq_processing/final_outs/", fastq)
  cellbender_seurat_dir <- list.files(path = cellbender_seurat_path, pattern = "_cellbender_seurat.rds", full.names = TRUE)
  cellbender_seurat <- readRDS(cellbender_seurat_dir)
  
  cells_in_sample <- Cells(cellbender_seurat)
  cells_to_keep_in_sample <- intersect(cells_in_sample, AD_cells)
  
  cellbender_seurat_filtered <- subset(cellbender_seurat, cells = cells_to_keep_in_sample)
  
  filtered_save_dir <- paste0("/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ad/", fastq, "_cellbender_seurat_filtered.rds")
  saveRDS(cellbender_seurat_filtered, file = filtered_save_dir)
}

########################################

# Gerrits FTD

FTD_fastq_list <- unique(sub(".*_", "", rownames(FTD_meta_filtered)))

# for testing: fastq <- "C4A"

for (fastq in FTD_fastq_list){
  cellbender_seurat_path <- paste0("/data/ADRD/gerrits_ftd_snRNA/fastq_processing/final_outs/", fastq)
  cellbender_seurat_dir <- list.files(path = cellbender_seurat_path, pattern = "_cellbender_seurat.rds", full.names = TRUE)
  cellbender_seurat <- readRDS(cellbender_seurat_dir)
  
  colnames(cellbender_seurat) <- paste0(colnames(cellbender_seurat), "_", fastq)
  
  cells_in_sample <- Cells(cellbender_seurat)
  cells_to_keep_in_sample <- intersect(cells_in_sample, FTD_cells)
  
  cellbender_seurat_filtered <- subset(cellbender_seurat, cells = cells_to_keep_in_sample)
  
  cellbender_seurat_filtered$sample <- fastq
  
  filtered_save_dir <- paste0("/data/ADRD/glia_across_NDDs/combined_data/post_qc_filter_individual_objects/ftd/", fastq, "_cellbender_seurat_filtered.rds")
  saveRDS(cellbender_seurat_filtered, file = filtered_save_dir)
}

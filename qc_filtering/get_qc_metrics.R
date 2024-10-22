library(Seurat)
library(tidyverse)
library(scCustomize)

########################################

# PD

fastq_list <- read.csv(file = "/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt", sep = "\t", header = F)
fastq_list <- fastq_list$V1

# loop through all the cellranger_seurat RDS files, pull relevant metadata

combined_metadata <- data.frame()

for (fastq in fastq_list){
  cellranger_seurat_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/", fastq, "/", fastq, "_cellranger_seurat.rds")
  cellranger_seurat <- readRDS(cellranger_seurat_dir)
  
  meta <- cellranger_seurat@meta.data[, c("donor", "batch", "region", "nCount_RNA", "nFeature_RNA", "pct.mito", "pct.ribo")]
  
  combined_metadata <- rbind(combined_metadata, meta)
}

# filter out cells from samples that had the same donor twice in the same pool

t <- as.data.frame(table(combined_metadata$batch))

combined_metadata_filtered <- combined_metadata %>%
  filter(!(batch == "PM-PD_Set4_1" & donor == "PM-MS_70978")) %>%
  filter(!(batch == "PM-PD_Set4_2" & donor == "PM-MS_70978")) %>%
  filter(!(batch == "PM-PD_Set5_1" & donor == "PM-MS_46475")) %>%
  filter(!(batch == "PM-PD_Set5_2" & donor == "PM-MS_46475")) %>%
  filter(!(batch == "PM-PD_Set7a_1" & donor == "PM-MS_55245")) %>%
  filter(!(batch == "PM-PD_Set7a_2" & donor == "PM-MS_55245"))

t2 <- as.data.frame(table(combined_metadata_filtered$batch))

saveRDS(combined_metadata_filtered, file = "/data/ADRD/glia_across_NDDs/amp_pd/combined_metadata_prefilter.rds")

########################################

# ALS/FTD

fastq_list <- read.csv(file = "/data/ADRD/ALSFTD_multiregion/fastq_processing/cleaned_SRRs.txt", sep = "\t", header = F)
fastq_list <- fastq_list$V1

# loop through objects & pull metadata

combined_metadata <- data.frame()

for (fastq in fastq_list){
  cellranger_seurat_path <- paste0("/data/ADRD/ALSFTD_multiregion/fastq_processing/final_outs/", fastq)
  cellranger_seurat_dir <- list.files(path = cellranger_seurat_path, pattern = "_cellranger_seurat.rds", full.names = TRUE)
  cellranger_seurat <- readRDS(cellranger_seurat_dir)
  
  meta <- cellranger_seurat@meta.data[, c("sample", "region", "nCount_RNA", "nFeature_RNA", "pct.mito", "pct.ribo")]
  
  combined_metadata <- rbind(combined_metadata, meta)
}

saveRDS(combined_metadata, file = "/data/ADRD/glia_across_NDDs/alsftd/combined_metadata_prefilter.rds")

########################################

# AD

fastq_list <- read.csv(file = "/data/ADRD/AD_multiregion/fastq_processing/round_1_fastqs.txt", sep = "\t", header = F)
fastq_list <- fastq_list$V1

# loop through objects & pull metadata

combined_metadata <- data.frame()

for (fastq in fastq_list){
  cellranger_seurat_path <- paste0("/data/ADRD/AD_multiregion/fastq_processing/final_outs/", fastq)
  cellranger_seurat_dir <- list.files(path = cellranger_seurat_path, pattern = "_cellranger_seurat.rds", full.names = TRUE)
  cellranger_seurat <- readRDS(cellranger_seurat_dir)
  
  meta <- cellranger_seurat@meta.data[, c("sample", "region", "nCount_RNA", "nFeature_RNA", "pct.mito", "pct.ribo")]
  
  combined_metadata <- rbind(combined_metadata, meta)
}

saveRDS(combined_metadata, file = "/data/ADRD/glia_across_NDDs/ad/combined_metadata_prefilter.rds")

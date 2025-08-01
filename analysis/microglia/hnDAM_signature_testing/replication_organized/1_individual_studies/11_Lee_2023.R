library(edgeR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(SeuratDisk)
library(EnvStats)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# Processed data from Lee et al., 2023, Sci Adv., downloaded from GEO (GSE148434)

# read in / merge sample objects

samples <- list.files("./analysis/microglia/DAM_signature_replication/individual_studies/Lee_2023/processed_data_from_study/", recursive = F, full.names = F)
# sample = samples[1]

seurat_list <- list()


for (sample in samples){
  sample <- gsub("_snRNA.molecule_info.h5", "", sample)
  print(paste0("Processing ", sample))
  cts <- Read10X_h5_GEO(paste0("./analysis/microglia/DAM_signature_replication/individual_studies/Lee_2023/processed_data_from_study/", sample, "_snRNA.molecule_info.h5"))
  seurat <- CreateSeuratObject(counts = cts)
  
  
  sample <- sub("^[^_]*_", "", sample)
  seurat$sample <- sample
  colnames(seurat) <- paste0(colnames(seurat), "_", sample)
  
  seurat_list[[sample]] <- seurat
}






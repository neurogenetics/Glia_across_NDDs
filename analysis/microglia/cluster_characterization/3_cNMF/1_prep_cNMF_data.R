library(tidyverse)
library(Seurat)
library(scCustomize)

setwd("/data/ADRD/glia_across_NDDs/")

set.seed(12345)

########################################

APPPD_cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/AMP-PD_microglia_sc_cts.rds")
Gerrits_cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Gerrits_microglia_sc_cts.rds")
Mathys_cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Mathys_microglia_sc_cts.rds")
Pineda_cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Pineda_microglia_sc_cts.rds")

celllevel_meta <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")

t <- celllevel_meta %>%
  group_by(dataset) %>%
  summarise(count = n())

########################################
########################################
########################################

# randomly sample 25000 cells per dataset for speed purposes

downsampled_cells <- celllevel_meta %>%
  rownames_to_column(var = "barcode") %>%
  group_by(dataset) %>%
  slice_sample(n = 25000, replace = FALSE) %>%
  ungroup() %>%
  column_to_rownames(var = "barcode")


AMPPD_cells_keep <- rownames(downsampled_cells)[downsampled_cells$dataset == "AMP-PD"]
Gerrits_cells_keep <- rownames(downsampled_cells)[downsampled_cells$dataset == "Gerrits"]
Mathys_cells_keep <- rownames(downsampled_cells)[downsampled_cells$dataset == "Mathys"]
Pineda_cells_keep <- rownames(downsampled_cells)[downsampled_cells$dataset == "Pineda"]


AMPPD_meta_ds <- downsampled_cells[downsampled_cells$dataset == "AMP-PD",]
Gerrits_meta_ds <- downsampled_cells[downsampled_cells$dataset == "Gerrits",]
Mathys_meta_ds <- downsampled_cells[downsampled_cells$dataset == "Mathys",]
Pineda_meta_ds <- downsampled_cells[downsampled_cells$dataset == "Pineda",]


AMPPD_meta_ds <- AMPPD_meta_ds[c("nCount_RNA", "nFeature_RNA", "donor", "region", "dataset", "pct.mito", "cluster_anno")]
Gerrits_meta_ds <- Gerrits_meta_ds[c("nCount_RNA", "nFeature_RNA", "donor", "region", "dataset", "pct.mito", "cluster_anno")]
Mathys_meta_ds <- Mathys_meta_ds[c("nCount_RNA", "nFeature_RNA", "donor", "region", "dataset", "pct.mito", "cluster_anno")]
Pineda_meta_ds <- Pineda_meta_ds[c("nCount_RNA", "nFeature_RNA", "donor", "region", "dataset", "pct.mito", "cluster_anno")]

########################################

# filter counts matrices for downsampled cells

APPPD_cts_ds <- APPPD_cts[,colnames(APPPD_cts) %in% AMPPD_cells_keep]
Gerrits_cts_ds <- Gerrits_cts[,colnames(Gerrits_cts) %in% Gerrits_cells_keep]
Mathys_cts_ds <- Mathys_cts[,colnames(Mathys_cts) %in% Mathys_cells_keep]
Pineda_cts_ds <- Pineda_cts[,colnames(Pineda_cts) %in% Pineda_cells_keep]

########################################

# make Seurat objects

seurat_AMPPD_cNMF <- CreateSeuratObject(counts = APPPD_cts_ds, meta.data = AMPPD_meta_ds)
seurat_Gerrits_cNMF <- CreateSeuratObject(counts = Gerrits_cts_ds, meta.data = Gerrits_meta_ds)
seurat_Mathys_cNMF <- CreateSeuratObject(counts = Mathys_cts_ds, meta.data = Mathys_meta_ds)
seurat_Pineda_cNMF <- CreateSeuratObject(counts = Pineda_cts_ds, meta.data = Pineda_meta_ds)


# save to .h5s

as.anndata(seurat_AMPPD_cNMF, 
           file_path = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/",
           file_name = "microglia_cNMF_AMPPD",
           main_layer = "counts", 
           other_layers = NULL)

as.anndata(seurat_Gerrits_cNMF, 
           file_path = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/",
           file_name = "microglia_cNMF_Gerrits",
           main_layer = "counts", 
           other_layers = NULL)

as.anndata(seurat_Mathys_cNMF, 
           file_path = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/",
           file_name = "microglia_cNMF_Mathys",
           main_layer = "counts", 
           other_layers = NULL)

as.anndata(seurat_Pineda_cNMF, 
           file_path = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/",
           file_name = "microglia_cNMF_Pineda",
           main_layer = "counts", 
           other_layers = NULL)

library(Seurat)
library(scCustomize)
library(tidyverse)
library(BPCells)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs")

########################################

PD_astrocytes <- readRDS("./combined_data/processed_celltypes_by_dataset/PD_astrocytes_allregions_filtered.rds")
AD_astrocytes <- readRDS("./combined_data/processed_celltypes_by_dataset/AD_astrocytes_allregions_filtered.rds")
ALS_astrocytes <- readRDS("./combined_data/processed_celltypes_by_dataset/ALS_astrocytes_allregions_filtered.rds")
FTD_astrocytes <- readRDS("./combined_data/processed_celltypes_by_dataset/FTD_astrocytes_allregions_filtered.rds")

########################################

# unify metadata across datasets

metadata_to_keep <- c("donor", "region", "pct.mito", "Donor")

PD_astrocytes@meta.data <- PD_astrocytes@meta.data[, colnames(PD_astrocytes@meta.data) %in% metadata_to_keep]
PD_astrocytes$dataset <- "AMP-PD"

ALS_astrocytes@meta.data <- ALS_astrocytes@meta.data[, colnames(ALS_astrocytes@meta.data) %in% metadata_to_keep]
ALS_astrocytes$dataset <- "Pineda"

AD_astrocytes@meta.data <- AD_astrocytes@meta.data[, colnames(AD_astrocytes@meta.data) %in% metadata_to_keep]
AD_astrocytes$dataset <- "Mathys"

FTD_astrocytes@meta.data <- FTD_astrocytes@meta.data[, colnames(FTD_astrocytes@meta.data) %in% metadata_to_keep]
FTD_astrocytes$dataset <- "Gerrits"
FTD_astrocytes@meta.data <- FTD_astrocytes@meta.data %>%
  dplyr::rename(donor = Donor)

########################################

# get raw counts matrices and filter genes expressed in 1% of all cells (or not -- testing both)

PD_counts <- GetAssayData(PD_astrocytes, assay = "RNA", layer = "counts")
AD_counts <- GetAssayData(AD_astrocytes, assay = "RNA", layer = "counts")
ALS_counts <- GetAssayData(ALS_astrocytes, assay = "RNA", layer = "counts")
FTD_counts <- GetAssayData(FTD_astrocytes, assay = "RNA", layer = "counts")

counts_combined <- cbind(PD_counts, AD_counts, ALS_counts, FTD_counts)
cell_count_threshold <- round(ncol(counts_combined) * 0.01)
counts_filtered <- counts_combined[rowSums(counts_combined > 0) >= cell_count_threshold, ]


PD_meta <- PD_astrocytes@meta.data
AD_meta <- AD_astrocytes@meta.data
ALS_meta <- ALS_astrocytes@meta.data
FTD_meta <- FTD_astrocytes@meta.data
meta_merged <- rbind(PD_meta, AD_meta, ALS_meta, FTD_meta)


astro <- CreateSeuratObject(counts = counts_combined, meta.data = meta_merged)
astro@meta.data <- astro@meta.data[, !colnames(astro@meta.data) %in% "orig.ident"]

astro_ctsfilter <- CreateSeuratObject(counts = counts_filtered, meta.data = meta_merged)
astro_ctsfilter@meta.data <- astro_ctsfilter@meta.data[, !colnames(astro_ctsfilter@meta.data) %in% "orig.ident"]

########################################
########################################
########################################

# write and read BPCells objects (for objects w/ filtered or unfiltered counts matrices)

write_matrix_dir(mat = astro[["RNA"]]$counts, 
                 dir = "./combined_data/testing_cts_filtering/astro_NO_ctsfilter_BPCells")
astro_mat <- open_matrix_dir(dir = "./combined_data/testing_cts_filtering/astro_NO_ctsfilter_BPCells")
astro_seurat <- CreateSeuratObject(counts = astro_mat, meta.data = astro@meta.data)
saveRDS(object = astro_seurat, 
        file = "./combined_data/testing_cts_filtering/astro_NO_ctsfilter_BPCells_seurat_RAW.rds")



write_matrix_dir(mat = astro_ctsfilter[["RNA"]]$counts, 
                 dir = "./combined_data/testing_cts_filtering/astro_YES_ctsfilter_BPCells")
astro_ctsfilter_mat <- open_matrix_dir(dir = "./combined_data/testing_cts_filtering/astro_YES_ctsfilter_BPCells")
astro_ctsfilter_seurat <- CreateSeuratObject(counts = astro_ctsfilter_mat, meta.data = astro_ctsfilter@meta.data)
saveRDS(object = astro_ctsfilter_seurat, 
        file = "./combined_data/testing_cts_filtering/astro_YES_ctsfilter_BPCells_seurat_RAW.rds")

########################################
########################################
########################################

# sketch datasets

astro <- readRDS("./combined_data/testing_cts_filtering/astro_NO_ctsfilter_BPCells_seurat_RAW.rds")
astro$dataset_region <- paste0(astro$dataset, "_", astro$region)
astro <- NormalizeData(astro)
astro[["RNA"]] <- split(astro[["RNA"]], f = astro$dataset_region)
astro <- FindVariableFeatures(astro)
astro_sketched <- SketchData(object = astro, ncells = 13000, method = "LeverageScore", sketched.assay = "sketch")
saveRDS(astro_sketched, file = "./combined_data/testing_cts_filtering/astro_sketched_13k_NO_ctsfilter_unprocessed.rds")



astro_ctsfilter <- readRDS("./combined_data/testing_cts_filtering/astro_YES_ctsfilter_BPCells_seurat_RAW.rds")
astro_ctsfilter$dataset_region <- paste0(astro_ctsfilter$dataset, "_", astro_ctsfilter$region)
astro_ctsfilter <- NormalizeData(astro_ctsfilter)
astro_ctsfilter[["RNA"]] <- split(astro_ctsfilter[["RNA"]], f = astro_ctsfilter$dataset_region)
astro_ctsfilter <- FindVariableFeatures(astro_ctsfilter)
astro_ctsfilter_sketched <- SketchData(object = astro_ctsfilter, ncells = 13000, method = "LeverageScore", sketched.assay = "sketch")
saveRDS(astro_ctsfilter_sketched, file = "./combined_data/testing_cts_filtering/astro_sketched_13k_YES_ctsfilter_unprocessed.rds")

########################################
########################################
########################################

set.seed(12345)

astro_sketched <- readRDS("./combined_data/testing_cts_filtering/astro_sketched_13k_NO_ctsfilter_unprocessed.rds")
DefaultAssay(astro_sketched) <- "sketch"
astro_sketched <- FindVariableFeatures(astro_sketched)
astro_sketched <- ScaleData(astro_sketched, vars.to.regress = "pct.mito")
astro_sketched <- RunPCA(astro_sketched)
astro_sketched[["sketch"]] <- JoinLayers(astro_sketched[["sketch"]])
astro_sketched <- RunHarmony(object = astro_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"), 
                             reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(astro_sketched, ndims = 50, reduction = "harmony")
astro_sketched <- FindNeighbors(astro_sketched, reduction = "harmony", dims = 1:25)
astro_sketched <- RunUMAP(astro_sketched, reduction = "harmony", dims = 1:25, return.model = T)
astro_sketched <- FindClusters(astro_sketched, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5))
astro_cluster_meta <- astro_sketched@meta.data
saveRDS(astro_cluster_meta, "./combined_data/testing_cts_filtering/metadata_for_clustree/astro_sketched_metadata_for_clustree_NO_ctsfilter.rds")
saveRDS(astro_sketched, file = "./combined_data/testing_cts_filtering/astro_sketched_clustered_NOT_projected_NO_ctsfilter.rds")



astro_ctsfilter_sketched <- readRDS("./combined_data/testing_cts_filtering/astro_sketched_13k_YES_ctsfilter_unprocessed.rds")
DefaultAssay(astro_ctsfilter_sketched) <- "sketch"
astro_ctsfilter_sketched <- FindVariableFeatures(astro_ctsfilter_sketched)
astro_ctsfilter_sketched <- ScaleData(astro_ctsfilter_sketched, vars.to.regress = "pct.mito")
astro_ctsfilter_sketched <- RunPCA(astro_ctsfilter_sketched)
astro_ctsfilter_sketched[["sketch"]] <- JoinLayers(astro_ctsfilter_sketched[["sketch"]])
astro_ctsfilter_sketched <- RunHarmony(object = astro_ctsfilter_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"), 
                                       reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(astro_ctsfilter_sketched, ndims = 50, reduction = "harmony")
astro_ctsfilter_sketched <- FindNeighbors(astro_ctsfilter_sketched, reduction = "harmony", dims = 1:25)
astro_ctsfilter_sketched <- RunUMAP(astro_ctsfilter_sketched, reduction = "harmony", dims = 1:25, return.model = T)
astro_ctsfilter_sketched <- FindClusters(astro_ctsfilter_sketched, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5))
astro_ctsfilter_cluster_meta <- astro_ctsfilter_sketched@meta.data
saveRDS(astro_ctsfilter_cluster_meta, "./combined_data/testing_cts_filtering/metadata_for_clustree/astro_sketched_metadata_for_clustree_YES_ctsfilter.rds")
saveRDS(astro_ctsfilter_sketched, file = "./combined_data/testing_cts_filtering/astro_sketched_clustered_NOT_projected_YES_ctsfilter.rds")

########################################
########################################
########################################

astro_sketched <- readRDS("./combined_data/testing_cts_filtering/astro_sketched_clustered_NOT_projected_NO_ctsfilter.rds")
astro_ctsfilter_sketched <- readRDS("./combined_data/testing_cts_filtering/astro_sketched_clustered_NOT_projected_YES_ctsfilter.rds")

res_list <- c("0.025", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5")


astro_markers_list <- list()
astro_avgexp_list <- list()

astro_ctsfilter_markers_list <- list()
astro_ctsfilter_avgexp_list <- list()


for (res in res_list){
  print(paste0("Finding markers and avgexp for res", res))
  
  Idents(astro_sketched) <- paste0("sketch_snn_res.", res)
  markers <- FindAllMarkers(astro_sketched, logfc.threshold = 0.5, min.pct = 0.2, only.pos = T)
  avgexp <- as.data.frame(AverageExpression(astro_sketched, assay = "RNA", layer = "counts"))
  
  astro_markers_list[[res]] <- markers
  astro_avgexp_list[[res]] <- avgexp
  
  
  Idents(astro_ctsfilter_sketched) <- paste0("sketch_snn_res.", res)
  markers_ctsfilter <- FindAllMarkers(astro_ctsfilter_sketched, logfc.threshold = 0.5, min.pct = 0.2, only.pos = T)
  avgexp_ctsfilter <- as.data.frame(AverageExpression(astro_ctsfilter_sketched, assay = "RNA", layer = "counts"))
  
  astro_ctsfilter_markers_list[[res]] <- markers_ctsfilter
  astro_ctsfilter_avgexp_list[[res]] <- avgexp_ctsfilter
}


saveRDS(astro_markers_list, file = "./combined_data/testing_cts_filtering/markers_avgexp/astro_markers_across_resolutions_NO_ctsfilter.rds")
saveRDS(astro_avgexp_list, file = "./combined_data/testing_cts_filtering/markers_avgexp/astro_avgexp_across_resolutions_NO_ctsfilter.rds")
saveRDS(astro_ctsfilter_markers_list, file = "./combined_data/testing_cts_filtering/markers_avgexp/astro_markers_across_resolutions_YES_ctsfilter.rds")
saveRDS(astro_ctsfilter_avgexp_list, file = "./combined_data/testing_cts_filtering/markers_avgexp/astro_avgexp_across_resolutions_YES_ctsfilter.rds")


astro_varfeats <- VariableFeatures(astro_sketched)
astro_ctsfilter_varfeats <- VariableFeatures(astro_ctsfilter_sketched)
saveRDS(astro_varfeats, file = "./combined_data/testing_cts_filtering/markers_avgexp/astro_varfeats_NO_ctsfilter.rds")
saveRDS(astro_ctsfilter_varfeats, file = "./combined_data/testing_cts_filtering/markers_avgexp/astro_varfeats_YES_ctsfilter.rds")


########################################
########################################
########################################

# going with no counts filter, res 0.15

set.seed(12345)
options(future.globals.maxSize = 8000 * 1024^2)

astro_sketched <- readRDS("./combined_data/testing_cts_filtering/astro_sketched_clustered_NOT_projected_NO_ctsfilter.rds")

astro_sketched[["sketch"]] <- split(astro_sketched[["sketch"]], f = astro_sketched$dataset_region)

astro_sketched <- ProjectIntegration(object = astro_sketched, sketched.assay = "sketch", 
                                     assay = "RNA", reduction = "harmony")

astro_sketched <- ProjectData(object = astro_sketched, sketched.assay = "sketch", assay = "RNA", 
                              sketched.reduction = "harmony.full",
                              full.reduction = "harmony.full", dims = 1:25, 
                              refdata = list(celltype.full = "sketch_snn_res.0.15"))

astro_sketched <- RunUMAP(astro_sketched, reduction = "harmony.full", dims = 1:25, 
                          reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(astro_sketched) <- "RNA"
Idents(astro_sketched) <- "celltype.full"

saveRDS(astro_sketched, file = "./combined_data/final_objects/astro_sketched_clustered_projected_NO_ctsfilter_res_0.15.rds")

########################################

astro_meta <- astro_sketched@meta.data
saveRDS(astro_meta, file = "./combined_data/final_objects/astro_projected_celllevel_meta.rds")

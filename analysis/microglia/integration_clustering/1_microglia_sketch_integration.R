library(Seurat)
library(scCustomize)
library(tidyverse)
library(BPCells)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs")

########################################

PD_microglia <- readRDS("./combined_data/processed_celltypes_by_dataset/PD_microglia_allregions_filtered.rds")
AD_microglia <- readRDS("./combined_data/processed_celltypes_by_dataset/AD_microglia_allregions_filtered.rds")
ALS_microglia <- readRDS("./combined_data/processed_celltypes_by_dataset/ALS_microglia_allregions_filtered.rds")
FTD_microglia <- readRDS("./combined_data/processed_celltypes_by_dataset/FTD_microglia_allregions_filtered.rds")

########################################

# unify metadata across datasets

metadata_to_keep <- c("donor", "region", "pct.mito", "Donor")

PD_microglia@meta.data <- PD_microglia@meta.data[, colnames(PD_microglia@meta.data) %in% metadata_to_keep]
PD_microglia$dataset <- "AMP-PD"

ALS_microglia@meta.data <- ALS_microglia@meta.data[, colnames(ALS_microglia@meta.data) %in% metadata_to_keep]
ALS_microglia$dataset <- "Pineda"

AD_microglia@meta.data <- AD_microglia@meta.data[, colnames(AD_microglia@meta.data) %in% metadata_to_keep]
AD_microglia$dataset <- "Mathys"

FTD_microglia@meta.data <- FTD_microglia@meta.data[, colnames(FTD_microglia@meta.data) %in% metadata_to_keep]
FTD_microglia$dataset <- "Gerrits"
FTD_microglia@meta.data <- FTD_microglia@meta.data %>%
  dplyr::rename(donor = Donor)

########################################

# get raw counts matrices and filter genes expressed in 1% of all cells (or not -- testing both)

PD_counts <- GetAssayData(PD_microglia, assay = "RNA", layer = "counts")
AD_counts <- GetAssayData(AD_microglia, assay = "RNA", layer = "counts")
ALS_counts <- GetAssayData(ALS_microglia, assay = "RNA", layer = "counts")
FTD_counts <- GetAssayData(FTD_microglia, assay = "RNA", layer = "counts")

counts_combined <- cbind(PD_counts, AD_counts, ALS_counts, FTD_counts)
cell_count_threshold <- round(ncol(counts_combined) * 0.01)
counts_filtered <- counts_combined[rowSums(counts_combined > 0) >= cell_count_threshold, ]


PD_meta <- PD_microglia@meta.data
AD_meta <- AD_microglia@meta.data
ALS_meta <- ALS_microglia@meta.data
FTD_meta <- FTD_microglia@meta.data
meta_merged <- rbind(PD_meta, AD_meta, ALS_meta, FTD_meta)


micro <- CreateSeuratObject(counts = counts_combined, meta.data = meta_merged)
micro@meta.data <- micro@meta.data[, !colnames(micro@meta.data) %in% "orig.ident"]

micro_ctsfilter <- CreateSeuratObject(counts = counts_filtered, meta.data = meta_merged)
micro_ctsfilter@meta.data <- micro_ctsfilter@meta.data[, !colnames(micro_ctsfilter@meta.data) %in% "orig.ident"]

########################################
########################################
########################################

# write and read BPCells objects (for objects w/ filtered or unfiltered counts matrices)

write_matrix_dir(mat = micro[["RNA"]]$counts, 
                 dir = "./combined_data/testing_cts_filtering/micro_NO_ctsfilter_BPCells")
micro_mat <- open_matrix_dir(dir = "./combined_data/testing_cts_filtering/micro_NO_ctsfilter_BPCells")
micro_seurat <- CreateSeuratObject(counts = micro_mat, meta.data = micro@meta.data)
saveRDS(object = micro_seurat, 
        file = "./combined_data/testing_cts_filtering/micro_NO_ctsfilter_BPCells_seurat_RAW.rds")



write_matrix_dir(mat = micro_ctsfilter[["RNA"]]$counts, 
                 dir = "./combined_data/testing_cts_filtering/micro_YES_ctsfilter_BPCells")
micro_ctsfilter_mat <- open_matrix_dir(dir = "./combined_data/testing_cts_filtering/micro_YES_ctsfilter_BPCells")
micro_ctsfilter_seurat <- CreateSeuratObject(counts = micro_ctsfilter_mat, meta.data = micro_ctsfilter@meta.data)
saveRDS(object = micro_ctsfilter_seurat, 
        file = "./combined_data/testing_cts_filtering/micro_YES_ctsfilter_BPCells_seurat_RAW.rds")

########################################
########################################
########################################

# sketch datasets

micro <- readRDS("./combined_data/testing_cts_filtering/micro_NO_ctsfilter_BPCells_seurat_RAW.rds")
micro$dataset_region <- paste0(micro$dataset, "_", micro$region)
micro <- NormalizeData(micro)
micro[["RNA"]] <- split(micro[["RNA"]], f = micro$dataset_region)
micro <- FindVariableFeatures(micro)
micro_sketched <- SketchData(object = micro, ncells = 6000, method = "LeverageScore", sketched.assay = "sketch")
saveRDS(micro_sketched, file = "./combined_data/testing_cts_filtering/micro_sketched_6k_NO_ctsfilter_unprocessed.rds")



micro_ctsfilter <- readRDS("./combined_data/testing_cts_filtering/micro_YES_ctsfilter_BPCells_seurat_RAW.rds")
micro_ctsfilter$dataset_region <- paste0(micro_ctsfilter$dataset, "_", micro_ctsfilter$region)
micro_ctsfilter <- NormalizeData(micro_ctsfilter)
micro_ctsfilter[["RNA"]] <- split(micro_ctsfilter[["RNA"]], f = micro_ctsfilter$dataset_region)
micro_ctsfilter <- FindVariableFeatures(micro_ctsfilter)
micro_ctsfilter_sketched <- SketchData(object = micro_ctsfilter, ncells = 6000, method = "LeverageScore", sketched.assay = "sketch")
saveRDS(micro_ctsfilter_sketched, file = "./combined_data/testing_cts_filtering/micro_sketched_6k_YES_ctsfilter_unprocessed.rds")

########################################
########################################
########################################

set.seed(12345)

micro_sketched <- readRDS("./combined_data/testing_cts_filtering/micro_sketched_6k_NO_ctsfilter_unprocessed.rds")
DefaultAssay(micro_sketched) <- "sketch"
micro_sketched <- FindVariableFeatures(micro_sketched)
micro_sketched <- ScaleData(micro_sketched, vars.to.regress = "pct.mito")
micro_sketched <- RunPCA(micro_sketched)
micro_sketched[["sketch"]] <- JoinLayers(micro_sketched[["sketch"]])
micro_sketched <- RunHarmony(object = micro_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"), 
                             reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(micro_sketched, ndims = 50, reduction = "harmony")
micro_sketched <- FindNeighbors(micro_sketched, reduction = "harmony", dims = 1:25)
micro_sketched <- RunUMAP(micro_sketched, reduction = "harmony", dims = 1:25, return.model = T)
micro_sketched <- FindClusters(micro_sketched, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5))
micro_cluster_meta <- micro_sketched@meta.data
saveRDS(micro_cluster_meta, "./combined_data/testing_cts_filtering/metadata_for_clustree/micro_sketched_metadata_for_clustree_NO_ctsfilter.rds")
saveRDS(micro_sketched, file = "./combined_data/testing_cts_filtering/micro_sketched_clustered_NOT_projected_NO_ctsfilter.rds")



micro_ctsfilter_sketched <- readRDS("./combined_data/testing_cts_filtering/micro_sketched_6k_YES_ctsfilter_unprocessed.rds")
DefaultAssay(micro_ctsfilter_sketched) <- "sketch"
micro_ctsfilter_sketched <- FindVariableFeatures(micro_ctsfilter_sketched)
micro_ctsfilter_sketched <- ScaleData(micro_ctsfilter_sketched, vars.to.regress = "pct.mito")
micro_ctsfilter_sketched <- RunPCA(micro_ctsfilter_sketched)
micro_ctsfilter_sketched[["sketch"]] <- JoinLayers(micro_ctsfilter_sketched[["sketch"]])
micro_ctsfilter_sketched <- RunHarmony(object = micro_ctsfilter_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"), 
                             reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(micro_ctsfilter_sketched, ndims = 50, reduction = "harmony")
micro_ctsfilter_sketched <- FindNeighbors(micro_ctsfilter_sketched, reduction = "harmony", dims = 1:25)
micro_ctsfilter_sketched <- RunUMAP(micro_ctsfilter_sketched, reduction = "harmony", dims = 1:25, return.model = T)
micro_ctsfilter_sketched <- FindClusters(micro_ctsfilter_sketched, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5))
micro_ctsfilter_cluster_meta <- micro_ctsfilter_sketched@meta.data
saveRDS(micro_ctsfilter_cluster_meta, "./combined_data/testing_cts_filtering/metadata_for_clustree/micro_sketched_metadata_for_clustree_YES_ctsfilter.rds")
saveRDS(micro_ctsfilter_sketched, file = "./combined_data/testing_cts_filtering/micro_sketched_clustered_NOT_projected_YES_ctsfilter.rds")

########################################
########################################
########################################

micro_sketched <- readRDS("./combined_data/testing_cts_filtering/micro_sketched_clustered_NOT_projected_NO_ctsfilter.rds")
micro_ctsfilter_sketched <- readRDS("./combined_data/testing_cts_filtering/micro_sketched_clustered_NOT_projected_YES_ctsfilter.rds")

res_list <- c("0.025", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5")


micro_markers_list <- list()
micro_avgexp_list <- list()

micro_ctsfilter_markers_list <- list()
micro_ctsfilter_avgexp_list <- list()


for (res in res_list){
  print(paste0("Finding markers and avgexp for res", res))
  
  Idents(micro_sketched) <- paste0("sketch_snn_res.", res)
  markers <- FindAllMarkers(micro_sketched, logfc.threshold = 0.5, min.pct = 0.2, only.pos = T)
  avgexp <- as.data.frame(AverageExpression(micro_sketched, assay = "RNA", layer = "counts"))
  
  micro_markers_list[[res]] <- markers
  micro_avgexp_list[[res]] <- avgexp
  
  
  Idents(micro_ctsfilter_sketched) <- paste0("sketch_snn_res.", res)
  markers_ctsfilter <- FindAllMarkers(micro_ctsfilter_sketched, logfc.threshold = 0.5, min.pct = 0.2, only.pos = T)
  avgexp_ctsfilter <- as.data.frame(AverageExpression(micro_ctsfilter_sketched, assay = "RNA", layer = "counts"))
  
  micro_ctsfilter_markers_list[[res]] <- markers_ctsfilter
  micro_ctsfilter_avgexp_list[[res]] <- avgexp_ctsfilter
}


saveRDS(micro_markers_list, file = "./combined_data/testing_cts_filtering/markers_avgexp/micro_markers_across_resolutions_NO_ctsfilter.rds")
saveRDS(micro_avgexp_list, file = "./combined_data/testing_cts_filtering/markers_avgexp/micro_avgexp_across_resolutions_NO_ctsfilter.rds")
saveRDS(micro_ctsfilter_markers_list, file = "./combined_data/testing_cts_filtering/markers_avgexp/micro_markers_across_resolutions_YES_ctsfilter.rds")
saveRDS(micro_ctsfilter_avgexp_list, file = "./combined_data/testing_cts_filtering/markers_avgexp/micro_avgexp_across_resolutions_YES_ctsfilter.rds")


micro_varfeats <- VariableFeatures(micro_sketched)
micro_ctsfilter_varfeats <- VariableFeatures(micro_ctsfilter_sketched)
saveRDS(micro_varfeats, file = "./combined_data/testing_cts_filtering/markers_avgexp/micro_varfeats_NO_ctsfilter.rds")
saveRDS(micro_ctsfilter_varfeats, file = "./combined_data/testing_cts_filtering/markers_avgexp/micro_varfeats_YES_ctsfilter.rds")


########################################
########################################
########################################

# going with no counts filter, res 0.25

set.seed(12345)

micro_sketched <- readRDS("./combined_data/testing_cts_filtering/micro_sketched_clustered_NOT_projected_NO_ctsfilter.rds")

micro_sketched[["sketch"]] <- split(micro_sketched[["sketch"]], f = micro_sketched$dataset_region)

micro_sketched <- ProjectIntegration(object = micro_sketched, sketched.assay = "sketch", 
                                     assay = "RNA", reduction = "harmony")

micro_sketched <- ProjectData(object = micro_sketched, sketched.assay = "sketch", assay = "RNA", 
                              sketched.reduction = "harmony.full",
                              full.reduction = "harmony.full", dims = 1:25, 
                              refdata = list(celltype.full = "sketch_snn_res.0.25"))

micro_sketched <- RunUMAP(micro_sketched, reduction = "harmony.full", dims = 1:25, 
                          reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(micro_sketched) <- "RNA"
Idents(micro_sketched) <- "celltype.full"

saveRDS(micro_sketched, file = "./combined_data/final_objects/micro_sketched_clustered_projected_NO_ctsfilter_res_0.25.rds")

library(Seurat)
library(tidyverse)
library(harmony)

########################################

######
# PD #
######

# functions

preprocess_pd_step_1 <- function(seurat_obj){
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  seurat_obj$batch <- gsub("PM-PD_", "", seurat_obj$batch)
  seurat_obj$pct.mito <- PercentageFeatureSet(seurat_obj, pattern = "^MT")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "pct.mito"))
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunHarmony(object = seurat_obj, reduction = "pca", group.by.vars = c("donor", "batch"),
                           reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
  
  return(seurat_obj)
}

preprocess_pd_step_2 <- function(seurat_obj, ndims, res){
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:ndims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:ndims, reduction.name = "umap_harmony")
  seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = "harmony_clusters")
  
  return(seurat_obj)
}

######

# running by region -- preprocess, check clusters, subset micro/astro, save

DMV_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/DMV_merged.rds")
DMV_merged <- preprocess_pd_step_1(DMV_merged)
ElbowPlot(DMV_merged, ndims = 50, reduction = "harmony")
DMV_merged <- preprocess_pd_step_2(DMV_merged, ndims = 30, res = 0.1)
Idents(DMV_merged) <- "harmony_clusters"
DimPlot(DMV_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(DMV_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
DMV_merged$region <- "DMV"
saveRDS(DMV_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/pd/DMV_all_preprocessed.rds")
DMV_microglia <- subset(DMV_merged, idents = "3")
DMV_astrocytes <- subset(DMV_merged, idents = "2")
saveRDS(DMV_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/DMV_microglia_raw.rds")
saveRDS(DMV_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/DMV_astrocytes_raw.rds")
DMV_pcs <- as.data.frame(DMV_merged@reductions$harmony@stdev)
DMV_varfeats <- as.data.frame(VariableFeatures(DMV_merged))
saveRDS(DMV_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/DMV_pcs.rds")
saveRDS(DMV_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/DMV_varfeats.rds")


GPi_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/GPi_merged.rds")
GPi_merged <- preprocess_pd_step_1(GPi_merged)
ElbowPlot(GPi_merged, ndims = 50, reduction = "harmony")
GPi_merged <- preprocess_pd_step_2(GPi_merged, ndims = 30, res = 0.1)
Idents(GPi_merged) <- "harmony_clusters"
DimPlot(GPi_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(GPi_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
GPi_merged$region <- "GPi"
saveRDS(GPi_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/pd/GPi_all_preprocessed.rds")
GPi_microglia <- subset(GPi_merged, idents = "2")
GPi_astrocytes <- subset(GPi_merged, idents = "4")
saveRDS(GPi_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/GPi_microglia_raw.rds")
saveRDS(GPi_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/GPi_astrocytes_raw.rds")
GPi_pcs <- as.data.frame(GPi_merged@reductions$harmony@stdev)
GPi_varfeats <- as.data.frame(VariableFeatures(GPi_merged))
saveRDS(GPi_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/GPi_pcs.rds")
saveRDS(GPi_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/GPi_varfeats.rds")


M1_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/M1_merged.rds")
M1_merged <- preprocess_pd_step_1(M1_merged)
ElbowPlot(M1_merged, ndims = 50, reduction = "harmony")
M1_merged <- preprocess_pd_step_2(M1_merged, ndims = 30, res = 0.1)
Idents(M1_merged) <- "harmony_clusters"
DimPlot(M1_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(M1_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
M1_merged$region <- "M1"
saveRDS(M1_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/pd/M1_all_preprocessed.rds")
M1_microglia <- subset(M1_merged, idents = "5")
M1_astrocytes <- subset(M1_merged, idents = c("2", "15"))
saveRDS(M1_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/M1_microglia_raw.rds")
saveRDS(M1_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/M1_astrocytes_raw.rds")
M1_pcs <- as.data.frame(M1_merged@reductions$harmony@stdev)
M1_varfeats <- as.data.frame(VariableFeatures(M1_merged))
saveRDS(M1_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/M1_pcs.rds")
saveRDS(M1_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/M1_varfeats.rds")


PFC_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/PFC_merged.rds")
PFC_merged <- preprocess_pd_step_1(PFC_merged)
ElbowPlot(PFC_merged, ndims = 50, reduction = "harmony")
PFC_merged <- preprocess_pd_step_2(PFC_merged, ndims = 30, res = 0.1)
Idents(PFC_merged) <- "harmony_clusters"
DimPlot(PFC_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(PFC_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
PFC_merged$region <- "PFC"
saveRDS(PFC_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/pd/PFC_all_preprocessed.rds")
PFC_microglia <- subset(PFC_merged, idents = "7")
PFC_astrocytes <- subset(PFC_merged, idents = c("2", "17"))
saveRDS(PFC_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/PFC_microglia_raw.rds")
saveRDS(PFC_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/PFC_astrocytes_raw.rds")
PFC_pcs <- as.data.frame(PFC_merged@reductions$harmony@stdev)
PFC_varfeats <- as.data.frame(VariableFeatures(PFC_merged))
saveRDS(PFC_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/PFC_pcs.rds")
saveRDS(PFC_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/PFC_varfeats.rds")


V1_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/pd/V1_merged.rds")
V1_merged <- preprocess_pd_step_1(V1_merged)
ElbowPlot(V1_merged, ndims = 50, reduction = "harmony")
V1_merged <- preprocess_pd_step_2(V1_merged, ndims = 30, res = 0.1)
Idents(V1_merged) <- "harmony_clusters"
DimPlot(V1_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(V1_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
V1_merged$region <- "V1"
saveRDS(V1_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/pd/V1_all_preprocessed.rds")
V1_microglia <- subset(V1_merged, idents = "6")
V1_astrocytes <- subset(V1_merged, idents = c("2", "17"))
saveRDS(V1_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/V1_microglia_raw.rds")
saveRDS(V1_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/pd/V1_astrocytes_raw.rds")
V1_pcs <- as.data.frame(V1_merged@reductions$harmony@stdev)
V1_varfeats <- as.data.frame(VariableFeatures(V1_merged))
saveRDS(V1_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/V1_pcs.rds")
saveRDS(V1_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/pd/V1_varfeats.rds")


########################################
########################################
########################################

###########
# ALS/FTD #
###########

# need to add collection date to metadata for integration

als_sample_meta <- readRDS("/data/ADRD/glia_across_NDDs/metadata/cleaned/als_sample_meta.rds")
colnames(als_sample_meta) <- c("srr_id", "collection_date")


# functions

preprocess_als_step_1 <- function(seurat_obj){
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  seurat_obj$pct.mito <- PercentageFeatureSet(seurat_obj, pattern = "^MT")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "pct.mito"))
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunHarmony(object = seurat_obj, reduction = "pca", group.by.vars = c("donor", "srr_id", "collection_date"),
                           reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
  
  return(seurat_obj)
}

preprocess_als_step_2 <- function(seurat_obj, ndims, res){
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:ndims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:ndims, reduction.name = "umap_harmony")
  seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = "harmony_clusters")
  
  return(seurat_obj)
}

######

# running by region -- preprocess, check clusters, subset micro/astro, save

M1_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/als/M1_merged.rds")
M1_merged$barcode <- colnames(M1_merged)
M1_meta <- M1_merged@meta.data 
M1_meta_merged <- M1_meta %>%
  left_join(als_sample_meta, by = "srr_id") %>%
  dplyr::rename(donor = sample) %>%
  dplyr::rename(donor_region = sample_region)
rownames(M1_meta_merged) <- M1_meta_merged$barcode
identical(rownames(M1_meta_merged), rownames(M1_meta))
M1_merged@meta.data <- M1_meta_merged
M1_merged <- preprocess_als_step_1(M1_merged)
ElbowPlot(M1_merged, ndims = 50, reduction = "harmony")
M1_merged <- preprocess_als_step_2(M1_merged, ndims = 30, res = 0.1)
Idents(M1_merged) <- "harmony_clusters"
DimPlot(M1_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(M1_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
saveRDS(M1_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/als/M1_all_preprocessed.rds")
M1_microglia <- subset(M1_merged, idents = "5")
M1_astrocytes <- subset(M1_merged, idents = c("2", "17"))
saveRDS(M1_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/als/M1_microglia_raw.rds")
saveRDS(M1_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/als/M1_astrocytes_raw.rds")
M1_pcs <- as.data.frame(M1_merged@reductions$harmony@stdev)
M1_varfeats <- as.data.frame(VariableFeatures(M1_merged))
saveRDS(M1_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/als/M1_pcs.rds")
saveRDS(M1_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/als/M1_varfeats.rds")


PFC_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/als/PFC_merged.rds")
PFC_merged$barcode <- colnames(PFC_merged)
PFC_meta <- PFC_merged@meta.data 
PFC_meta_merged <- PFC_meta %>%
  left_join(als_sample_meta, by = "srr_id") %>%
  dplyr::rename(donor = sample) %>%
  dplyr::rename(donor_region = sample_region)
rownames(PFC_meta_merged) <- PFC_meta_merged$barcode
identical(rownames(PFC_meta_merged), rownames(PFC_meta))
PFC_merged@meta.data <- PFC_meta_merged
PFC_merged <- preprocess_als_step_1(PFC_merged)
ElbowPlot(PFC_merged, ndims = 50, reduction = "harmony")
PFC_merged <- preprocess_als_step_2(PFC_merged, ndims = 30, res = 0.1)
Idents(PFC_merged) <- "harmony_clusters"
DimPlot(PFC_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(PFC_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
saveRDS(PFC_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/als/PFC_all_preprocessed.rds")
PFC_microglia <- subset(PFC_merged, idents = "7")
PFC_astrocytes <- subset(PFC_merged, idents = "2")
saveRDS(PFC_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/als/PFC_microglia_raw.rds")
saveRDS(PFC_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/als/PFC_astrocytes_raw.rds")
PFC_pcs <- as.data.frame(PFC_merged@reductions$harmony@stdev)
PFC_varfeats <- as.data.frame(VariableFeatures(PFC_merged))
saveRDS(PFC_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/als/PFC_pcs.rds")
saveRDS(PFC_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/als/PFC_varfeats.rds")

########################################
########################################
########################################

######
# AD #
######

# functions

preprocess_ad_step_1 <- function(seurat_obj){
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  seurat_obj$pct.mito <- PercentageFeatureSet(seurat_obj, pattern = "^MT")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "pct.mito"))
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunHarmony(object = seurat_obj, reduction = "pca", group.by.vars = c("donor"),
                           reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
  
  return(seurat_obj)
}

preprocess_ad_step_2 <- function(seurat_obj, ndims, res){
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:ndims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:ndims, reduction.name = "umap_harmony")
  seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = "harmony_clusters")
  
  return(seurat_obj)
}

######

# running by region -- preprocess, check clusters, subset micro/astro, save

AnG_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/AnG_merged.rds")
AnG_merged@meta.data <- AnG_merged@meta.data %>%
  dplyr::rename(donor = sample) %>%
  dplyr::rename(donor_region = sample_region)
AnG_merged <- preprocess_ad_step_1(AnG_merged)
ElbowPlot(AnG_merged, ndims = 50, reduction = "harmony")
AnG_merged <- preprocess_ad_step_2(AnG_merged, ndims = 30, res = 0.1)
Idents(AnG_merged) <- "harmony_clusters"
DimPlot(AnG_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(AnG_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
saveRDS(AnG_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/ad/AnG_all_preprocessed.rds")
AnG_microglia <- subset(AnG_merged, idents = "6")
AnG_astrocytes <- subset(AnG_merged, idents = "2")
saveRDS(AnG_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/AnG_microglia_raw.rds")
saveRDS(AnG_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/AnG_astrocytes_raw.rds")
AnG_pcs <- as.data.frame(AnG_merged@reductions$harmony@stdev)
AnG_varfeats <- as.data.frame(VariableFeatures(AnG_merged))
saveRDS(AnG_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/AnG_pcs.rds")
saveRDS(AnG_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/AnG_varfeats.rds")


HIP_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/HIP_merged.rds")
HIP_merged@meta.data <- HIP_merged@meta.data %>%
  dplyr::rename(donor = sample) %>%
  dplyr::rename(donor_region = sample_region)
HIP_merged <- preprocess_ad_step_1(HIP_merged)
ElbowPlot(HIP_merged, ndims = 50, reduction = "harmony")
HIP_merged <- preprocess_ad_step_2(HIP_merged, ndims = 30, res = 0.1)
Idents(HIP_merged) <- "harmony_clusters"
DimPlot(HIP_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(HIP_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
saveRDS(HIP_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/ad/HIP_all_preprocessed.rds")
HIP_microglia <- subset(HIP_merged, idents = "2")
HIP_astrocytes <- subset(HIP_merged, idents = c("1", "14"))
saveRDS(HIP_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/HIP_microglia_raw.rds")
saveRDS(HIP_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/HIP_astrocytes_raw.rds")
HIP_pcs <- as.data.frame(HIP_merged@reductions$harmony@stdev)
HIP_varfeats <- as.data.frame(VariableFeatures(HIP_merged))
saveRDS(HIP_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/HIP_pcs.rds")
saveRDS(HIP_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/HIP_varfeats.rds")


EC_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/EC_merged.rds")
EC_merged@meta.data <- EC_merged@meta.data %>%
  dplyr::rename(donor = sample) %>%
  dplyr::rename(donor_region = sample_region)
EC_merged <- preprocess_ad_step_1(EC_merged)
ElbowPlot(EC_merged, ndims = 50, reduction = "harmony")
EC_merged <- preprocess_ad_step_2(EC_merged, ndims = 30, res = 0.1)
Idents(EC_merged) <- "harmony_clusters"
DimPlot(EC_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(EC_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
saveRDS(EC_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/ad/EC_all_preprocessed.rds")
EC_microglia <- subset(EC_merged, idents = "2")
EC_astrocytes <- subset(EC_merged, idents = c("1", "12"))
saveRDS(EC_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/EC_microglia_raw.rds")
saveRDS(EC_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/EC_astrocytes_raw.rds")
EC_pcs <- as.data.frame(EC_merged@reductions$harmony@stdev)
EC_varfeats <- as.data.frame(VariableFeatures(EC_merged))
saveRDS(EC_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/EC_pcs.rds")
saveRDS(EC_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/EC_varfeats.rds")


PFC_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/PFC_merged.rds")
PFC_merged@meta.data <- PFC_merged@meta.data %>%
  dplyr::rename(donor = sample) %>%
  dplyr::rename(donor_region = sample_region)
PFC_merged <- preprocess_ad_step_1(PFC_merged)
ElbowPlot(PFC_merged, ndims = 50, reduction = "harmony")
PFC_merged <- preprocess_ad_step_2(PFC_merged, ndims = 30, res = 0.1)
Idents(PFC_merged) <- "harmony_clusters"
DimPlot(PFC_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(PFC_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
saveRDS(PFC_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/ad/PFC_all_preprocessed.rds")
PFC_microglia <- subset(PFC_merged, idents = "8")
PFC_astrocytes <- subset(PFC_merged, idents = "2")
saveRDS(PFC_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/PFC_microglia_raw.rds")
saveRDS(PFC_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/PFC_astrocytes_raw.rds")
PFC_pcs <- as.data.frame(PFC_merged@reductions$harmony@stdev)
PFC_varfeats <- as.data.frame(VariableFeatures(PFC_merged))
saveRDS(PFC_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/PFC_pcs.rds")
saveRDS(PFC_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/PFC_varfeats.rds")


TH_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/TH_merged.rds")
TH_merged@meta.data <- TH_merged@meta.data %>%
  dplyr::rename(donor = sample) %>%
  dplyr::rename(donor_region = sample_region)
TH_merged <- preprocess_ad_step_1(TH_merged)
ElbowPlot(TH_merged, ndims = 50, reduction = "harmony")
TH_merged <- preprocess_ad_step_2(TH_merged, ndims = 30, res = 0.1)
Idents(TH_merged) <- "harmony_clusters"
DimPlot(TH_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(TH_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
saveRDS(TH_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/ad/TH_all_preprocessed.rds")
TH_microglia <- subset(TH_merged, idents = "3")
TH_astrocytes <- subset(TH_merged, idents = c("2", "6"))
saveRDS(TH_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/TH_microglia_raw.rds")
saveRDS(TH_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/TH_astrocytes_raw.rds")
TH_pcs <- as.data.frame(TH_merged@reductions$harmony@stdev)
TH_varfeats <- as.data.frame(VariableFeatures(TH_merged))
saveRDS(TH_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/TH_pcs.rds")
saveRDS(TH_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/TH_varfeats.rds")


MTG_merged <- readRDS("/data/ADRD/glia_across_NDDs/combined_data/merged_regions/ad/MTG_merged.rds")
MTG_merged@meta.data <- MTG_merged@meta.data %>%
  dplyr::rename(donor = sample) %>%
  dplyr::rename(donor_region = sample_region)
MTG_merged <- preprocess_ad_step_1(MTG_merged)
ElbowPlot(MTG_merged, ndims = 50, reduction = "harmony")
MTG_merged <- preprocess_ad_step_2(MTG_merged, ndims = 30, res = 0.1)
Idents(MTG_merged) <- "harmony_clusters"
DimPlot(MTG_merged, reduction = "umap_harmony", label = T) 
FeaturePlot(MTG_merged, features = c("RBFOX3", "SLC17A7", "GAD1", "P2RY12", "AQP4", "MOG"), ncol = 3, label = T)
saveRDS(MTG_merged, "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/ad/MTG_all_preprocessed.rds")
MTG_microglia <- subset(MTG_merged, idents = "5")
MTG_astrocytes <- subset(MTG_merged, idents = "2")
saveRDS(MTG_microglia, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/MTG_microglia_raw.rds")
saveRDS(MTG_astrocytes, "/data/ADRD/glia_across_NDDs/combined_data/subset_celltypes_by_region/ad/MTG_astrocytes_raw.rds")
MTG_pcs <- as.data.frame(MTG_merged@reductions$harmony@stdev)
MTG_varfeats <- as.data.frame(VariableFeatures(MTG_merged))
saveRDS(MTG_pcs, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/MTG_pcs.rds")
saveRDS(MTG_varfeats, file = "/data/ADRD/glia_across_NDDs/combined_data/processed_regions/varfeats_umapdims_from_processing/ad/MTG_varfeats.rds")

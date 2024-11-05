library(Seurat)
library(tidyverse)
library(scCustomize)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs/combined_data")

########################################

# PD

##################

# functions for processing

preprocess_pd_subsets_1 <- function(seurat_obj){
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "pct.mito"))
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunHarmony(object = seurat_obj, reduction = "pca", group.by.vars = c("donor", "batch"), 
                           reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
  
  return(seurat_obj)
}

preprocess_pd_subsets_2 <- function(seurat_obj, ndims, res, cluster_name){
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:ndims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:ndims, reduction.name = "umap_harmony")
  seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = cluster_name)
  
  return(seurat_obj)
}

##################

# microglia

PD_DMV_micro <- readRDS("./subset_celltypes_by_region/pd/DMV_microglia_raw.rds")
PD_GPi_micro <- readRDS("./subset_celltypes_by_region/pd/GPi_microglia_raw.rds")
PD_M1_micro <- readRDS("./subset_celltypes_by_region/pd/M1_microglia_raw.rds")
PD_PFC_micro <- readRDS("./subset_celltypes_by_region/pd/PFC_microglia_raw.rds")
PD_V1_micro <- readRDS("./subset_celltypes_by_region/pd/V1_microglia_raw.rds")

PD_micro_list <- c(PD_DMV_micro, PD_GPi_micro, PD_M1_micro, PD_PFC_micro, PD_V1_micro)
PD_micro_merged <- Merge_Seurat_List(PD_micro_list, merge.data = T)
PD_micro_merged <- preprocess_pd_subsets_1(PD_micro_merged)
ElbowPlot(PD_micro_merged, ndims = 50, reduction = "harmony")
PD_micro_merged <- preprocess_pd_subsets_2(PD_micro_merged, ndims = 20, res = 0.1, cluster_name = "harmony_microglia_clusters")

# check clusters and remove any that are not microglia

Idents(PD_micro_merged) <- "harmony_microglia_clusters"
DimPlot(PD_micro_merged, reduction = "umap_harmony", label = T, split.by = "region") 
FeaturePlot(PD_micro_merged, features = c("P2RY12", "MRC1", "MOG", "CD3E"), label = T)
PD_micro_markers <- FindAllMarkers(PD_micro_merged, logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
# 3 expresses oligodendrocyte marker MOG
# 8 are T-cells (CD3+)

PD_micro_filtered <- subset(PD_micro_merged, idents = c("0", "1", "2", "4", "5", "6", "7"))

# filter out donor x region w/ <50 cells

PD_micro_donor_counts <- as.data.frame(table(PD_micro_filtered$donor, PD_micro_filtered$region))
PD_micro_donor_counts <- PD_micro_donor_counts %>%
  arrange(Freq)
PD_micro_donor_counts$donor_region <- paste0(PD_micro_donor_counts$Var1, "_", PD_micro_donor_counts$Var2)
PD_donor_regions_to_keep <- PD_micro_donor_counts$donor_region[92:nrow(PD_micro_donor_counts)]

PD_micro_filtered$donor_region <- paste0(PD_micro_filtered$donor, "_", PD_micro_filtered$region)
PD_micro_filtered <- subset(PD_micro_filtered, subset = donor_region %in% PD_donor_regions_to_keep)

saveRDS(PD_micro_markers, file = "./processed_celltypes_by_dataset/markers_from_filtering/PD_micro_markers_prefilter.rds")
saveRDS(PD_micro_filtered, file = "./processed_celltypes_by_dataset/PD_microglia_allregions_filtered.rds")

##################

# astrocytes

PD_DMV_astro <- readRDS("./subset_celltypes_by_region/pd/DMV_astrocytes_raw.rds")
PD_GPi_astro <- readRDS("./subset_celltypes_by_region/pd/GPi_astrocytes_raw.rds")
PD_M1_astro <- readRDS("./subset_celltypes_by_region/pd/M1_astrocytes_raw.rds")
PD_PFC_astro <- readRDS("./subset_celltypes_by_region/pd/PFC_astrocytes_raw.rds")
PD_V1_astro <- readRDS("./subset_celltypes_by_region/pd/V1_astrocytes_raw.rds")

PD_astro_list <- c(PD_DMV_astro, PD_GPi_astro, PD_M1_astro, PD_PFC_astro, PD_V1_astro)
PD_astro_merged <- Merge_Seurat_List(PD_astro_list, merge.data = T)
PD_astro_merged <- preprocess_pd_subsets_1(PD_astro_merged)
ElbowPlot(PD_astro_merged, ndims = 50, reduction = "harmony")
PD_astro_merged <- preprocess_pd_subsets_2(PD_astro_merged, ndims = 30, res = 0.1, cluster_name = "harmony_astrocytes_clusters")

# check clusters and remove any that are not astrocytes

Idents(PD_astro_merged) <- "harmony_astrocytes_clusters"
DimPlot(PD_astro_merged, reduction = "umap_harmony", label = T, split.by = "region") 
FeaturePlot(PD_astro_merged, features = c("AQP4", "GFAP", "MOG", "P2RY12"), label = T)
PD_astro_markers <- FindAllMarkers(PD_astro_merged, logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
# 4 expresses oligodendrocyte markers (MOG, IL1RAPL1, ST18)
# 5 express microglial markers (DOCK8, CSF1R)
# 7 express OPC markers (KCNIP4, RBFOX1)

PD_astro_filtered <- subset(PD_astro_merged, idents = c("0", "1", "2", "3", "6", "8", "9"))

# filter out donor x region w/ <50 cells

PD_astro_donor_counts <- as.data.frame(table(PD_astro_filtered$donor, PD_astro_filtered$region))
PD_astro_donor_counts <- PD_astro_donor_counts %>%
  arrange(Freq)
PD_astro_donor_counts$donor_region <- paste0(PD_astro_donor_counts$Var1, "_", PD_astro_donor_counts$Var2)
PD_donor_regions_to_keep <- PD_astro_donor_counts$donor_region[43:nrow(PD_astro_donor_counts)]

PD_astro_filtered$donor_region <- paste0(PD_astro_filtered$donor, "_", PD_astro_filtered$region)
PD_astro_filtered <- subset(PD_astro_filtered, subset = donor_region %in% PD_donor_regions_to_keep)

saveRDS(PD_astro_markers, file = "./processed_celltypes_by_dataset/markers_from_filtering/PD_astro_markers_prefilter.rds")
saveRDS(PD_astro_filtered, file = "./processed_celltypes_by_dataset/PD_astrocytes_allregions_filtered.rds")

########################################
########################################
########################################

# ALS

##################

# functions for processing

preprocess_als_subsets_1 <- function(seurat_obj){
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "pct.mito"))
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunHarmony(object = seurat_obj, reduction = "pca", group.by.vars = c("donor", "srr_id", "collection_date"),
                           reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
  
  return(seurat_obj)
}

preprocess_als_subsets_2 <- function(seurat_obj, ndims, res, cluster_name){
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:ndims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:ndims, reduction.name = "umap_harmony")
  seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = cluster_name)
  
  return(seurat_obj)
}

##################

# microglia

ALS_M1_micro <- readRDS("./subset_celltypes_by_region/als/M1_microglia_raw.rds")
ALS_PFC_micro <- readRDS("./subset_celltypes_by_region/als/PFC_microglia_raw.rds")

ALS_micro_list <- c(ALS_M1_micro, ALS_PFC_micro)
ALS_micro_merged <- Merge_Seurat_List(ALS_micro_list, merge.data = T)
ALS_micro_merged <- preprocess_als_subsets_1(ALS_micro_merged)
ElbowPlot(ALS_micro_merged, ndims = 50, reduction = "harmony")
ALS_micro_merged <- preprocess_als_subsets_2(ALS_micro_merged, ndims = 20, res = 0.1, cluster_name = "harmony_microglia_clusters")

# check clusters and remove any that are not microglia

Idents(ALS_micro_merged) <- "harmony_microglia_clusters"
DimPlot(ALS_micro_merged, reduction = "umap_harmony", label = T, split.by = "region") 
FeaturePlot(ALS_micro_merged, features = c("P2RY12", "MRC1", "MOG", "CD3E"), label = T)
ALS_micro_markers <- FindAllMarkers(ALS_micro_merged, logfc.threshold = 0.5, min.pct = 0.3, only.pos = T)
# 3 expresses oligodendrocyte marker MOG
# 8 are T-cells (CD3+)

ALS_micro_filtered <- subset(ALS_micro_merged, idents = c("0", "1", "4", "5", "6"))

# filter out donor x region w/ <50 cells

ALS_micro_donor_counts <- as.data.frame(table(ALS_micro_filtered$donor, ALS_micro_filtered$region))
ALS_micro_donor_counts <- ALS_micro_donor_counts %>%
  arrange(Freq)
ALS_micro_donor_counts$donor_region <- paste0(ALS_micro_donor_counts$Var1, "_", ALS_micro_donor_counts$Var2)
ALS_donor_regions_to_keep <- ALS_micro_donor_counts$donor_region[22:nrow(ALS_micro_donor_counts)]

ALS_micro_filtered$donor_region <- paste0(ALS_micro_filtered$donor, "_", ALS_micro_filtered$region)
ALS_micro_filtered <- subset(ALS_micro_filtered, subset = donor_region %in% ALS_donor_regions_to_keep)

saveRDS(ALS_micro_markers, file = "./processed_celltypes_by_dataset/markers_from_filtering/ALS_micro_markers_prefilter.rds")
saveRDS(ALS_micro_filtered, file = "./processed_celltypes_by_dataset/ALS_microglia_allregions_filtered.rds")

##################

# astrocytes

ALS_M1_astro <- readRDS("./subset_celltypes_by_region/als/M1_astrocytes_raw.rds")
ALS_PFC_astro <- readRDS("./subset_celltypes_by_region/als/PFC_astrocytes_raw.rds")

ALS_astro_list <- c(ALS_M1_astro, ALS_PFC_astro)
ALS_astro_merged <- Merge_Seurat_List(ALS_astro_list, merge.data = T)
ALS_astro_merged <- preprocess_als_subsets_1(ALS_astro_merged)
ElbowPlot(ALS_astro_merged, ndims = 50, reduction = "harmony")
ALS_astro_merged <- preprocess_als_subsets_2(ALS_astro_merged, ndims = 20, res = 0.1, cluster_name = "harmony_astrocytes_clusters")

# check clusters and remove any that are not astrocytes

Idents(ALS_astro_merged) <- "harmony_astrocytes_clusters"
DimPlot(ALS_astro_merged, reduction = "umap_harmony", label = T, split.by = "region") 
FeaturePlot(ALS_astro_merged, features = c("AQP4", "GFAP", "MOG", "P2RY12"), label = T)
ALS_astro_markers <- FindAllMarkers(ALS_astro_merged, logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
# 4 expresses microglial markers
# 5 express oligo markers
# 6 very high in mito genes

ALS_astro_filtered <- subset(ALS_astro_merged, idents = c("0", "1", "2", "3"))

# filter out donor x region w/ <50 cells

ALS_astro_donor_counts <- as.data.frame(table(ALS_astro_filtered$donor, ALS_astro_filtered$region))
ALS_astro_donor_counts <- ALS_astro_donor_counts %>%
  arrange(Freq)
ALS_astro_donor_counts$donor_region <- paste0(ALS_astro_donor_counts$Var1, "_", ALS_astro_donor_counts$Var2)
ALS_donor_regions_to_keep <- ALS_astro_donor_counts$donor_region[6:nrow(ALS_astro_donor_counts)]

ALS_astro_filtered$donor_region <- paste0(ALS_astro_filtered$donor, "_", ALS_astro_filtered$region)
ALS_astro_filtered <- subset(ALS_astro_filtered, subset = donor_region %in% ALS_donor_regions_to_keep)

saveRDS(ALS_astro_markers, file = "./processed_celltypes_by_dataset/markers_from_filtering/ALS_astro_markers_prefilter.rds")
saveRDS(ALS_astro_filtered, file = "./processed_celltypes_by_dataset/ALS_astrocytes_allregions_filtered.rds")


########################################
########################################
########################################

# AD

##################

# functions for processing

preprocess_ad_subsets_1 <- function(seurat_obj){
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "pct.mito"))
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunHarmony(object = seurat_obj, reduction = "pca", group.by.vars = c("donor"), 
                           reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
  
  return(seurat_obj)
}

preprocess_ad_subsets_2 <- function(seurat_obj, ndims, res, cluster_name){
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:ndims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:ndims, reduction.name = "umap_harmony")
  seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = cluster_name)
  
  return(seurat_obj)
}

##################

# microglia

AD_TH_micro <- readRDS("./subset_celltypes_by_region/ad/TH_microglia_raw.rds")
AD_EC_micro <- readRDS("./subset_celltypes_by_region/ad/EC_microglia_raw.rds")
AD_HIP_micro <- readRDS("./subset_celltypes_by_region/ad/HIP_microglia_raw.rds")
AD_PFC_micro <- readRDS("./subset_celltypes_by_region/ad/PFC_microglia_raw.rds")
AD_MTG_micro <- readRDS("./subset_celltypes_by_region/ad/MTG_microglia_raw.rds")
AD_AnG_micro <- readRDS("./subset_celltypes_by_region/ad/AnG_microglia_raw.rds")

AD_micro_list <- c(AD_TH_micro, AD_EC_micro, AD_HIP_micro, AD_PFC_micro, AD_MTG_micro, AD_AnG_micro)
AD_micro_merged <- Merge_Seurat_List(AD_micro_list, merge.data = T)
AD_micro_merged <- preprocess_ad_subsets_1(AD_micro_merged)
ElbowPlot(AD_micro_merged, ndims = 50, reduction = "harmony")
AD_micro_merged <- preprocess_ad_subsets_2(AD_micro_merged, ndims = 20, res = 0.1, cluster_name = "harmony_microglia_clusters")

# check clusters and remove any that are not microglia

Idents(AD_micro_merged) <- "harmony_microglia_clusters"
DimPlot(AD_micro_merged, reduction = "umap_harmony", label = T, split.by = "region") 
FeaturePlot(AD_micro_merged, features = c("P2RY12", "MRC1", "MOG", "CD3E"), label = T)
AD_micro_markers <- FindAllMarkers(AD_micro_merged, logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
# 2, 7 expresses oligodendrocyte markers 
# 4 are T-cells/peripheral immune

AD_micro_filtered <- subset(AD_micro_merged, idents = c("0", "1", "6"))

# filter out donor x region w/ <50 cells

AD_micro_donor_counts <- as.data.frame(table(AD_micro_filtered$donor, AD_micro_filtered$region))
AD_micro_donor_counts <- AD_micro_donor_counts %>%
  arrange(Freq)
AD_micro_donor_counts$donor_region <- paste0(AD_micro_donor_counts$Var1, "_", AD_micro_donor_counts$Var2)
AD_donor_regions_to_keep <- AD_micro_donor_counts$donor_region[25:nrow(AD_micro_donor_counts)]

AD_micro_filtered$donor_region <- paste0(AD_micro_filtered$donor, "_", AD_micro_filtered$region)
AD_micro_filtered <- subset(AD_micro_filtered, subset = donor_region %in% AD_donor_regions_to_keep)

saveRDS(AD_micro_markers, file = "./processed_celltypes_by_dataset/markers_from_filtering/AD_micro_markers_prefilter.rds")
saveRDS(AD_micro_filtered, file = "./processed_celltypes_by_dataset/AD_microglia_allregions_filtered.rds")

##################

# astrocytes

AD_TH_astro <- readRDS("./subset_celltypes_by_region/ad/TH_astrocytes_raw.rds")
AD_EC_astro <- readRDS("./subset_celltypes_by_region/ad/EC_astrocytes_raw.rds")
AD_HIP_astro <- readRDS("./subset_celltypes_by_region/ad/HIP_astrocytes_raw.rds")
AD_PFC_astro <- readRDS("./subset_celltypes_by_region/ad/PFC_astrocytes_raw.rds")
AD_MTG_astro <- readRDS("./subset_celltypes_by_region/ad/MTG_astrocytes_raw.rds")
AD_AnG_astro <- readRDS("./subset_celltypes_by_region/ad/AnG_astrocytes_raw.rds")

AD_astro_list <- c(AD_TH_astro, AD_EC_astro, AD_HIP_astro, AD_PFC_astro, AD_MTG_astro, AD_AnG_astro)
AD_astro_merged <- Merge_Seurat_List(AD_astro_list, merge.data = T)
AD_astro_merged <- preprocess_ad_subsets_1(AD_astro_merged)
ElbowPlot(AD_astro_merged, ndims = 50, reduction = "harmony")
AD_astro_merged <- preprocess_ad_subsets_2(AD_astro_merged, ndims = 30, res = 0.1, cluster_name = "harmony_astrocytes_clusters")

# check clusters and remove any that are not astrocytes

Idents(AD_astro_merged) <- "harmony_astrocytes_clusters"
DimPlot(AD_astro_merged, reduction = "umap_harmony", label = T, split.by = "region") 
FeaturePlot(AD_astro_merged, features = c("AQP4", "GFAP", "MOG", "P2RY12"), label = T)
AD_astro_markers <- FindAllMarkers(AD_astro_merged, logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
# 3 are ependymal cells (SPAG17)
# 4 express oligo markers
# 6 express microglial markers (DOCK8)

AD_astro_filtered <- subset(AD_astro_merged, idents = c("0", "1", "2", "5"))

# filter out donor x region w/ <50 cells

AD_astro_donor_counts <- as.data.frame(table(AD_astro_filtered$donor, AD_astro_filtered$region))
AD_astro_donor_counts <- AD_astro_donor_counts %>%
  arrange(Freq)
AD_astro_donor_counts$donor_region <- paste0(AD_astro_donor_counts$Var1, "_", AD_astro_donor_counts$Var2)
AD_donor_regions_to_keep <- AD_astro_donor_counts$donor_region[16:nrow(AD_astro_donor_counts)]

AD_astro_filtered$donor_region <- paste0(AD_astro_filtered$donor, "_", AD_astro_filtered$region)
AD_astro_filtered <- subset(AD_astro_filtered, subset = donor_region %in% AD_donor_regions_to_keep)

saveRDS(AD_astro_markers, file = "./processed_celltypes_by_dataset/markers_from_filtering/AD_astro_markers_prefilter.rds")
saveRDS(AD_astro_filtered, file = "./processed_celltypes_by_dataset/AD_astrocytes_allregions_filtered.rds")

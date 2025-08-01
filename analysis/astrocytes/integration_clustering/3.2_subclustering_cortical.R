library(Seurat)
library(scCustomize)
library(tidyverse)
library(harmony)
library(clustree)

setwd("/data/ADRD/glia_across_NDDs")

########################################

all_astros_meta <- readRDS("./analysis/astrocytes/all_astrocytes_celllevel_metadata.rds")
AMPPD_astro_cts <- readRDS("./analysis/astrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_astrocytes_sc_cts.rds")
Mathys_astro_cts <- readRDS("./analysis/astrocytes/differential_expression/single_cell_cts_matrices/Mathys_astrocytes_sc_cts.rds")
Gerrits_astro_cts <- readRDS("./analysis/astrocytes/differential_expression/single_cell_cts_matrices/Gerrits_astrocytes_sc_cts.rds")
Pineda_astro_cts <- readRDS("./analysis/astrocytes/differential_expression/single_cell_cts_matrices/Pineda_astrocytes_sc_cts.rds")

cortical_astros_meta <- all_astros_meta[all_astros_meta$region %in% c("FC/PFC", "M1", "OC/V1", "EC", "HIP", "TC/MTG", "AnG"), ]
cortical_astros_meta <- cortical_astros_meta[c("nCount_RNA", "nFeature_RNA", "donor", "region", "pct.mito", "dataset")]

AMPPD_astro_cts <- AMPPD_astro_cts[, colnames(AMPPD_astro_cts) %in% rownames(cortical_astros_meta)]
Mathys_astro_cts <- Mathys_astro_cts[, colnames(Mathys_astro_cts) %in% rownames(cortical_astros_meta)]
Gerrits_astro_cts <- Gerrits_astro_cts[, colnames(Gerrits_astro_cts) %in% rownames(cortical_astros_meta)]
Pineda_astro_cts <- Pineda_astro_cts[, colnames(Pineda_astro_cts) %in% rownames(cortical_astros_meta)]

cts_merged <- cbind(AMPPD_astro_cts, Mathys_astro_cts, Gerrits_astro_cts, Pineda_astro_cts)


cortical_astros <- CreateSeuratObject(counts = cts_merged, meta.data = cortical_astros_meta)
cortical_astros[["RNA"]] <- JoinLayers(cortical_astros[["RNA"]])

########################################

# sketch data

cortical_astros$dataset_region <- paste0(cortical_astros$dataset, "_", cortical_astros$region)
cortical_astros[["RNA"]] <- split(cortical_astros[["RNA"]], f = cortical_astros$dataset_region)
cortical_astros <- NormalizeData(cortical_astros)
cortical_astros <- FindVariableFeatures(cortical_astros)
cortical_astros_sketched <- SketchData(object = cortical_astros, ncells = 13000, method = "LeverageScore", sketched.assay = "sketch")


DefaultAssay(cortical_astros_sketched) <- "sketch"
cortical_astros_sketched <- FindVariableFeatures(cortical_astros_sketched)
cortical_astros_sketched <- ScaleData(cortical_astros_sketched, vars.to.regress = "pct.mito")
cortical_astros_sketched <- RunPCA(cortical_astros_sketched)
cortical_astros_sketched[["sketch"]] <- JoinLayers(cortical_astros_sketched[["sketch"]])
cortical_astros_sketched <- RunHarmony(object = cortical_astros_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"),
                                          reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(cortical_astros_sketched, ndims = 50, reduction = "harmony")
cortical_astros_sketched <- FindNeighbors(cortical_astros_sketched, reduction = "harmony", dims = 1:15)
cortical_astros_sketched <- RunUMAP(cortical_astros_sketched, reduction = "harmony", dims = 1:15, return.model = T)
cortical_astros_sketched <- FindClusters(cortical_astros_sketched, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 
                                                                                        0.25, 0.3, 0.35, 0.4, 0.45, 0.5))

saveRDS(cortical_astros_sketched, file = "./analysis/astrocytes/seurat_objects/cortical_astrocytes_sketched_clustered_NOT_projected.rds")

########################################

DimPlot_scCustom(cortical_astros_sketched, group.by = "sketch_snn_res.0.2", split.by = "dataset")

DotPlot_scCustom(cortical_astros_sketched, 
                 features = c("AQP4", "GFAP", "PLP1", "MOG", "MBP", "SPP1"), group.by = "sketch_snn_res.0.25")

DotPlot_scCustom(cortical_astros_sketched, 
                 features = c("HSP90AA1", "HSPH1", "DNAJB1"), group.by = "sketch_snn_res.0.15")


Idents(cortical_astros_sketched) <- "sketch_snn_res.0.15"
Cluster_Highlight_Plot(cortical_astros_sketched, cluster_name = "4")

t <- cortical_astros_sketched@meta.data %>%
  group_by(sketch_snn_res.0.2) %>%
  summarise(umis = mean(nCount_RNA),
            mito = mean(pct.mito),
            count = n())

markers_0.15 <- FindAllMarkers(cortical_astros_sketched, group.by = "sketch_snn_res.0.15", only.pos = T, min.pct = 0.4, logfc.threshold = 0.5)
markers_0.2 <- FindAllMarkers(cortical_astros_sketched, group.by = "sketch_snn_res.0.2", only.pos = T, min.pct = 0.4, logfc.threshold = 0.5)
markers_0.25 <- FindAllMarkers(cortical_astros_sketched, group.by = "sketch_snn_res.0.25", only.pos = T, min.pct = 0.4, logfc.threshold = 0.5)

# c7 on res0.2 are ependymal cells (SPAG17)

FeaturePlot_scCustom(cortical_astros_sketched, features = "GRIA1")

########################################

cortical_astros_sketched <- readRDS("./analysis/astrocytes/seurat_objects/cortical_astrocytes_sketched_clustered_NOT_projected.rds")

# re-project at res0.2 to remove cluster 7

cortical_astros_sketched$dataset_region <- paste0(cortical_astros_sketched$dataset, "_", cortical_astros_sketched$region)
cortical_astros_sketched[["sketch"]] <- split(cortical_astros_sketched[["sketch"]], f = cortical_astros_sketched$dataset_region)

cortical_astros_sketched <- ProjectIntegration(object = cortical_astros_sketched, sketched.assay = "sketch",
                                                  assay = "RNA", reduction = "harmony")

options(future.globals.maxSize = 8000 * 1024^2)
cortical_astros_sketched <- ProjectData(object = cortical_astros_sketched, sketched.assay = "sketch", assay = "RNA",
                                           sketched.reduction = "harmony.full",
                                           full.reduction = "harmony.full", dims = 1:15,
                                           refdata = list(celltype.full = "sketch_snn_res.0.2"))

cortical_astros_sketched <- RunUMAP(cortical_astros_sketched, reduction = "harmony.full", dims = 1:15,
                                       reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(cortical_astros_sketched) <- "RNA"
Idents(cortical_astros_sketched) <- "celltype.full"

########################################

# removing cluster 7 

cells_to_keep <- colnames(cortical_astros_sketched)[!cortical_astros_sketched$celltype.full == "7"]
cortical_astros_sketched_filtered <- subset(cortical_astros_sketched, cells = cells_to_keep)

########################################
########################################
########################################

cortical_astros_sketched_filtered[["RNA"]] <- JoinLayers(cortical_astros_sketched_filtered[["RNA"]])

########################################

# sketch data

cortical_astros_sketched_filtered$dataset_region <- paste0(cortical_astros_sketched_filtered$dataset, "_", cortical_astros_sketched_filtered$region)
cortical_astros_sketched_filtered[["RNA"]] <- split(cortical_astros_sketched_filtered[["RNA"]], f = cortical_astros_sketched_filtered$dataset_region)
cortical_astros_sketched_filtered <- NormalizeData(cortical_astros_sketched_filtered)
cortical_astros_sketched_filtered <- FindVariableFeatures(cortical_astros_sketched_filtered)
seurat <- SketchData(object = cortical_astros_sketched_filtered, ncells = 13000, method = "LeverageScore", sketched.assay = "sketch")


DefaultAssay(seurat) <- "sketch"
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat, vars.to.regress = "pct.mito")
seurat <- RunPCA(seurat)
seurat[["sketch"]] <- JoinLayers(seurat[["sketch"]])
seurat <- RunHarmony(object = seurat, reduction = "pca", group.by.vars = c("dataset", "donor"),
                     reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(seurat, ndims = 50, reduction = "harmony")
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:15)
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:15, return.model = T)
seurat <- FindClusters(seurat, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5))

saveRDS(seurat, file = "./analysis/astrocytes/seurat_objects/cortical_astrocytes_sketched_filtered_clustered_NOT_projected.rds")

########################################

seurat <- readRDS("./analysis/astrocytes/seurat_objects/cortical_astrocytes_sketched_filtered_clustered_NOT_projected.rds")

seurat <- FindClusters(seurat, resolution = c(0.125))

########################################

DimPlot_scCustom(seurat, group.by = "sketch_snn_res.0.1")


markers_0.1 <- FindAllMarkers(seurat, group.by = "sketch_snn_res.0.1", min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
markers_0.125 <- FindAllMarkers(seurat, group.by = "sketch_snn_res.0.125", min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
markers_0.15 <- FindAllMarkers(seurat, group.by = "sketch_snn_res.0.15", min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)


FeaturePlot_scCustom(seurat, features = "CD44")

# going with res0.15

########################################

seurat$dataset_region <- paste0(seurat$dataset, "_", seurat$region)
seurat[["sketch"]] <- split(seurat[["sketch"]], f = seurat$dataset_region)

seurat <- ProjectIntegration(object = seurat, sketched.assay = "sketch",
                                               assay = "RNA", reduction = "harmony")

options(future.globals.maxSize = 8000 * 1024^2)
seurat <- ProjectData(object = seurat, sketched.assay = "sketch", assay = "RNA",
                                        sketched.reduction = "harmony.full",
                                        full.reduction = "harmony.full", dims = 1:15,
                                        refdata = list(celltype.full = "sketch_snn_res.0.15"))

seurat <- RunUMAP(seurat, reduction = "harmony.full", dims = 1:15,
                                    reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(seurat) <- "RNA"
Idents(seurat) <- "celltype.full"

saveRDS(seurat, file = "./analysis/astrocytes/seurat_objects/cortical_astrocytes_sketched_filtered_clustered_projected_0.15.rds")

library(Seurat)
library(scCustomize)
library(tidyverse)
library(harmony)
library(clustree)

setwd("/data/ADRD/glia_across_NDDs")

########################################

# get cells from subcortical regions (AMP-PD: DMV, GPi; Mathys: TH)

all_astros_meta <- readRDS("./analysis/astrocytes/all_astrocytes_celllevel_metadata.rds")
AMPPD_astro_cts <- readRDS("./analysis/astrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_astrocytes_sc_cts.rds")
Mathys_astro_cts <- readRDS("./analysis/astrocytes/differential_expression/single_cell_cts_matrices/Mathys_astrocytes_sc_cts.rds")

subcortical_astros_meta <- all_astros_meta[all_astros_meta$region %in% c("DMV", "GPi", "TH"), ]
subcortical_astros_meta <- subcortical_astros_meta[c("nCount_RNA", "nFeature_RNA", "donor", "region", "pct.mito", "dataset")]

AMPPD_astro_cts <- AMPPD_astro_cts[, colnames(AMPPD_astro_cts) %in% rownames(subcortical_astros_meta)]
Mathys_astro_cts <- Mathys_astro_cts[, colnames(Mathys_astro_cts) %in% rownames(subcortical_astros_meta)]

cts_merged <- cbind(AMPPD_astro_cts, Mathys_astro_cts)


subcortical_astros <- CreateSeuratObject(counts = cts_merged, meta.data = subcortical_astros_meta)
subcortical_astros[["RNA"]] <- JoinLayers(subcortical_astros[["RNA"]])

########################################
########################################
########################################

# sketch data

subcortical_astros$dataset_region <- paste0(subcortical_astros$dataset, "_", subcortical_astros$region)
subcortical_astros[["RNA"]] <- split(subcortical_astros[["RNA"]], f = subcortical_astros$dataset_region)
subcortical_astros <- NormalizeData(subcortical_astros)
subcortical_astros <- FindVariableFeatures(subcortical_astros)
subcortical_astros_sketched <- SketchData(object = subcortical_astros, ncells = 20000, method = "LeverageScore", sketched.assay = "sketch")


DefaultAssay(subcortical_astros_sketched) <- "sketch"
subcortical_astros_sketched <- FindVariableFeatures(subcortical_astros_sketched)
subcortical_astros_sketched <- ScaleData(subcortical_astros_sketched, vars.to.regress = "pct.mito")
subcortical_astros_sketched <- RunPCA(subcortical_astros_sketched)
subcortical_astros_sketched[["sketch"]] <- JoinLayers(subcortical_astros_sketched[["sketch"]])
subcortical_astros_sketched <- RunHarmony(object = subcortical_astros_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"),
                                          reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(subcortical_astros_sketched, ndims = 50, reduction = "harmony")
subcortical_astros_sketched <- FindNeighbors(subcortical_astros_sketched, reduction = "harmony", dims = 1:25)
subcortical_astros_sketched <- RunUMAP(subcortical_astros_sketched, reduction = "harmony", dims = 1:25, return.model = T)
subcortical_astros_sketched <- FindClusters(subcortical_astros_sketched, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 
                                                                                        0.25, 0.3, 0.35, 0.4, 0.45, 0.5))

saveRDS(subcortical_astros_sketched, file = "./analysis/astrocytes/seurat_objects/subcortical_astrocytes_sketched_clustered_NOT_projected.rds")

########################################
########################################
########################################

subcortical_astros_sketched <- readRDS("./analysis/astrocytes/seurat_objects/subcortical_astrocytes_sketched_clustered_NOT_projected.rds")

DimPlot_scCustom(subcortical_astros_sketched, group.by = "sketch_snn_res.0.2", split.by = "region")
DimPlot_scCustom(subcortical_astros_sketched, group.by = "sketch_snn_res.0.25", split.by = "region")

markers_0.2 <- FindAllMarkers(subcortical_astros_sketched, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T, group.by = "sketch_snn_res.0.2")
markers_0.25 <- FindAllMarkers(subcortical_astros_sketched, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T, group.by = "sketch_snn_res.0.25")

FeaturePlot_scCustom(subcortical_astros_sketched, features = "GRIN2A")
DotPlot_scCustom(subcortical_astros_sketched, features = c("AQP4", "GFAP", "PLP1", "MOG", "MBP", "SPP1"), group.by = "celltype.full")

t <- subcortical_astros_sketched@meta.data %>%
  group_by(sketch_snn_res.0.2) %>%
  summarise(avg_umis = mean(nCount_RNA), avg_mito = mean(pct.mito))

t2 <- subcortical_astros_sketched@meta.data %>%
  group_by(dataset, region, celltype.full) %>%
  summarise(counts = n())

########################################

# re-project at res0.2

subcortical_astros_sketched$dataset_region <- paste0(subcortical_astros_sketched$dataset, "_", subcortical_astros_sketched$region)
subcortical_astros_sketched[["sketch"]] <- split(subcortical_astros_sketched[["sketch"]], f = subcortical_astros_sketched$dataset_region)

subcortical_astros_sketched <- ProjectIntegration(object = subcortical_astros_sketched, sketched.assay = "sketch",
                                                  assay = "RNA", reduction = "harmony")

options(future.globals.maxSize = 8000 * 1024^2)
subcortical_astros_sketched <- ProjectData(object = subcortical_astros_sketched, sketched.assay = "sketch", assay = "RNA",
                                           sketched.reduction = "harmony.full",
                                           full.reduction = "harmony.full", dims = 1:25,
                                           refdata = list(celltype.full = "sketch_snn_res.0.2"))

subcortical_astros_sketched <- RunUMAP(subcortical_astros_sketched, reduction = "harmony.full", dims = 1:25,
                                       reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(subcortical_astros_sketched) <- "RNA"
Idents(subcortical_astros_sketched) <- "celltype.full"

DimPlot_scCustom(subcortical_astros_sketched)

########################################

# removing clusters 6 and 7 -- 6 = doublets, 7 = low UMIs and very high mitochondrial

cells_to_keep <- colnames(subcortical_astros_sketched)[!subcortical_astros_sketched$celltype.full %in% c("6", "7")]
seurat <- subset(subcortical_astros_sketched, cells = cells_to_keep)

########################################
########################################
########################################

seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])

########################################

# sketch data

seurat$dataset_region <- paste0(seurat$dataset, "_", seurat$region)
seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$dataset_region)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- SketchData(object = seurat, ncells = 20000, method = "LeverageScore", sketched.assay = "sketch")


DefaultAssay(seurat) <- "sketch"
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat, vars.to.regress = "pct.mito")
seurat <- RunPCA(seurat)
seurat[["sketch"]] <- JoinLayers(seurat[["sketch"]])
seurat <- RunHarmony(object = seurat, reduction = "pca", group.by.vars = c("dataset", "donor"),
                     reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(seurat, ndims = 50, reduction = "harmony")
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:25)
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:25, return.model = T)
seurat <- FindClusters(seurat, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5))

########################################

DimPlot_scCustom(seurat, group.by = "sketch_snn_res.0.2", split.by = "region")


markers_0.1 <- FindAllMarkers(seurat, group.by = "sketch_snn_res.0.1", logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
markers_0.15 <- FindAllMarkers(seurat, group.by = "sketch_snn_res.0.15", logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
markers_0.2 <- FindAllMarkers(seurat, group.by = "sketch_snn_res.0.2", logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)

meta <- na.omit(seurat@meta.data)
p <- clustree(meta, prefix = "sketch_snn_res.")


DotPlot_scCustom(seurat, group.by = "sketch_snn_res.0.2",
                 features = unique(c("GFAP", "AQP4",
                              "CABLES1", "SLC6A11", "FRMD4A",
                              "L3MBTL4", "CCDC85A", "LINC00609",
                              "NRXN3", "CTNNA3", "KSR2", 
                              "MDGA2", "SHISA6", "GRIA4", "ST3GAL6",
                              "DNAH6", "ARMC3", "LRRC9", 
                              "ARMC3", "CADPS", "NAMPT", "ZFYVE28", "CHI3L1",
                              "EPHA6", "GRM3", "ADGRV1", "GRIA2",
                              "SERPINH1", "HSPH1", "DNAJB1",
                              "SLIT2", "GRIA1", "KIRREL3"))) +
  RotatedAxis()

FeaturePlot_scCustom(seurat, features = "GFAP")

########################################

# re-project at res0.2

seurat$dataset_region <- paste0(seurat$dataset, "_", seurat$region)
seurat[["sketch"]] <- split(seurat[["sketch"]], f = seurat$dataset_region)

seurat <- ProjectIntegration(object = seurat, sketched.assay = "sketch",
                                                  assay = "RNA", reduction = "harmony")

options(future.globals.maxSize = 8000 * 1024^2)
seurat <- ProjectData(object = seurat, sketched.assay = "sketch", assay = "RNA",
                                           sketched.reduction = "harmony.full",
                                           full.reduction = "harmony.full", dims = 1:25,
                                           refdata = list(celltype.full = "sketch_snn_res.0.2"))

seurat <- RunUMAP(seurat, reduction = "harmony.full", dims = 1:25,
                                       reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(seurat) <- "RNA"
Idents(seurat) <- "celltype.full"

DimPlot_scCustom(seurat)

saveRDS(seurat, file = "./analysis/astrocytes/seurat_objects/subcortical_astrocytes_sketched_filtered_clustered_projected_0.2.rds")

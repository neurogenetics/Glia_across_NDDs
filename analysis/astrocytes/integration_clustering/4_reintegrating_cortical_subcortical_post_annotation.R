# library(Seurat) # manually load 5.2.1!!!
library(scCustomize)
library(tidyverse)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs/")

########################################

cortical <- readRDS("./analysis/astrocytes/seurat_objects/cortical_astrocytes_clustered_ANNOTATED.rds")
subcortical <- readRDS("./analysis/astrocytes/seurat_objects/subcortical_astrocytes_clustered_ANNOTATED.rds")

########################################

cortical_cts <- LayerData(cortical, layer = "counts", assay = "RNA")
subcortical_cts <- LayerData(subcortical, layer = "counts", assay = "RNA")

cts_merged <- cbind(cortical_cts, subcortical_cts)


cortical_meta <- cortical@meta.data %>%
  dplyr::select(nCount_RNA, nFeature_RNA, donor, region, pct.mito, dataset, dataset_region, cluster_anno)
subcortical_meta <- subcortical@meta.data %>%
  dplyr::select(nCount_RNA, nFeature_RNA, donor, region, pct.mito, dataset, dataset_region, cluster_anno)

meta_merged <- rbind(cortical_meta, subcortical_meta)


astros <- CreateSeuratObject(counts = cts_merged, meta.data = meta_merged)

########################################

# !!!! Need to make sure to use Seurat 5.2.1 for this!!! 5.3.0 is broken as of right now, throws an error on SketchData (4/28/2025)

astros[["RNA"]] <- split(astros[["RNA"]], f = astros$dataset_region)
astros <- NormalizeData(astros)
astros <- FindVariableFeatures(astros)
astros_sketched <- SketchData(object = astros, ncells = 10000, method = "LeverageScore", sketched.assay = "sketch")


DefaultAssay(astros_sketched) <- "sketch"
astros_sketched <- FindVariableFeatures(astros_sketched)
astros_sketched <- ScaleData(astros_sketched, vars.to.regress = "pct.mito")
astros_sketched <- RunPCA(astros_sketched)
astros_sketched[["sketch"]] <- JoinLayers(astros_sketched[["sketch"]])
astros_sketched <- RunHarmony(object = astros_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"),
                                          reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(astros_sketched, ndims = 50, reduction = "harmony")
astros_sketched <- FindNeighbors(astros_sketched, reduction = "harmony", dims = 1:20)
astros_sketched <- RunUMAP(astros_sketched, reduction = "harmony", dims = 1:20, return.model = T)
astros_sketched <- FindClusters(astros_sketched, resolution = c(0.1, 0.15, 0.2, 0.25, 0.3))

saveRDS(astros_sketched, file = "./analysis/astrocytes/seurat_objects/all_astros_sketched_clustered_post_annotation_NOT_projected.rds")

########################################

astros_sketched <- readRDS("./analysis/astrocytes/seurat_objects/all_astros_sketched_clustered_post_annotation_NOT_projected.rds")

# re-project at res0.15

astros_sketched$dataset_region <- paste0(astros_sketched$dataset, "_", astros_sketched$region)
astros_sketched[["sketch"]] <- split(astros_sketched[["sketch"]], f = astros_sketched$dataset_region)

astros_sketched <- ProjectIntegration(object = astros_sketched, sketched.assay = "sketch",
                                      assay = "RNA", reduction = "harmony")

options(future.globals.maxSize = 8000 * 1024^2)
astros_sketched <- ProjectData(object = astros_sketched, sketched.assay = "sketch", assay = "RNA",
                               sketched.reduction = "harmony.full", full.reduction = "harmony.full", dims = 1:20,
                               refdata = list(celltype.full = "sketch_snn_res.0.15"))

astros_sketched <- RunUMAP(astros_sketched, reduction = "harmony.full", dims = 1:20,
                           reduction.name = "umap.full", reduction.key = "UMAPfull_")

saveRDS(astros_sketched, file = "./analysis/astrocytes/seurat_objects/all_astros_merged_post_annotation_clustered_projected.rds")

library(Seurat)
library(scCustomize)
library(tidyverse)
library(harmony)
library(gdata)
library(clustree)

setwd("/data/ADRD/glia_across_NDDs")

########################################
########################################
########################################

PD_DMV_oligos_cts <- readRDS("./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_DMV_oligos_cts.rds")
PD_GPi_oligos_cts <- readRDS("./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_GPi_oligos_cts.rds")
PD_M1_oligos_cts <- readRDS("./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_M1_oligos_cts.rds")
PD_PFC_oligos_cts <- readRDS("./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_PFC_oligos_cts.rds")
PD_V1_oligos_cts <- readRDS("./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AMP-PD_V1_oligos_cts.rds")

AD_oligos_merged_cts <- readRDS("./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/AD_oligos_merged_cts.rds")

ALS_oligos_merged_cts <- readRDS("./analysis/oligodendrocytes/differential_expression/single_cell_cts_matrices/ALS_oligos_merged_cts.rds")

########################################

PD_DMV_oligos_meta <- readRDS("./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_DMV_oligos_meta.rds")
PD_GPi_oligos_meta <- readRDS("./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_GPi_oligos_meta.rds")
PD_M1_oligos_meta <- readRDS("./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_M1_oligos_meta.rds")
PD_PFC_oligos_meta <- readRDS("./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_PFC_oligos_meta.rds")
PD_V1_oligos_meta <- readRDS("./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AMP-PD_V1_oligos_meta.rds")

AD_oligos_merged_meta <- readRDS("./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/AD_oligos_merged_meta.rds")

ALS_oligos_merged_meta <- readRDS("./analysis/oligodendrocytes/seurat_objects/celllevel_meta_prior_to_merging/ALS_oligos_merged_meta.rds")

########################################

PD_DMV_oligos_meta <- PD_DMV_oligos_meta[c("nCount_RNA", "nFeature_RNA", "pct.mito", "donor", "region")]
PD_DMV_oligos_meta$dataset <- "AMP-PD"

PD_GPi_oligos_meta <- PD_GPi_oligos_meta[c("nCount_RNA", "nFeature_RNA", "pct.mito", "donor", "region")]
PD_GPi_oligos_meta$dataset <- "AMP-PD"

PD_M1_oligos_meta <- PD_M1_oligos_meta[c("nCount_RNA", "nFeature_RNA", "pct.mito", "donor", "region")]
PD_M1_oligos_meta$dataset <- "AMP-PD"

PD_PFC_oligos_meta <- PD_PFC_oligos_meta[c("nCount_RNA", "nFeature_RNA", "pct.mito", "donor", "region")]
PD_PFC_oligos_meta$dataset <- "AMP-PD"

PD_V1_oligos_meta <- PD_V1_oligos_meta[c("nCount_RNA", "nFeature_RNA", "pct.mito", "donor", "region")]
PD_V1_oligos_meta$dataset <- "AMP-PD"

AD_oligos_merged_meta <- AD_oligos_merged_meta[c("nCount_RNA", "nFeature_RNA", "pct.mito", "donor", "region")]
AD_oligos_merged_meta$dataset <- "Mathys"

ALS_oligos_merged_meta <- ALS_oligos_merged_meta[c("nCount_RNA", "nFeature_RNA", "pct.mito", "donor", "region")]
ALS_oligos_merged_meta$dataset <- "Pineda"

########################################
########################################
########################################

PD_DMV_oligos_seurat <- CreateSeuratObject(counts = PD_DMV_oligos_cts, meta.data = PD_DMV_oligos_meta)
PD_GPi_oligos_seurat <- CreateSeuratObject(counts = PD_GPi_oligos_cts, meta.data = PD_GPi_oligos_meta)
PD_M1_oligos_seurat <- CreateSeuratObject(counts = PD_M1_oligos_cts, meta.data = PD_M1_oligos_meta)
PD_PFC_oligos_seurat <- CreateSeuratObject(counts = PD_PFC_oligos_cts, meta.data = PD_PFC_oligos_meta)
PD_V1_oligos_seurat <- CreateSeuratObject(counts = PD_V1_oligos_cts, meta.data = PD_V1_oligos_meta)

AD_oligos_merged_seurat <- CreateSeuratObject(counts = AD_oligos_merged_cts, meta.data = AD_oligos_merged_meta)
AD_oligos_merged_seurat[["RNA"]] <- split(AD_oligos_merged_seurat[["RNA"]], f = AD_oligos_merged_seurat$region)

ALS_oligos_merged_seurat <- CreateSeuratObject(counts = ALS_oligos_merged_cts, meta.data = ALS_oligos_merged_meta)
ALS_oligos_merged_seurat[["RNA"]] <- split(ALS_oligos_merged_seurat[["RNA"]], f = ALS_oligos_merged_seurat$region)


rm(PD_DMV_oligos_cts)
rm(PD_GPi_oligos_cts)
rm(PD_M1_oligos_cts)
rm(PD_PFC_oligos_cts)
rm(PD_V1_oligos_cts)
rm(AD_oligos_merged_cts)
rm(ALS_oligos_merged_cts)

########################################

oligos_seurat_list <- list(PD_DMV_oligos_seurat, PD_GPi_oligos_seurat, PD_M1_oligos_seurat, PD_PFC_oligos_seurat,
                           PD_V1_oligos_seurat, AD_oligos_merged_seurat, ALS_oligos_merged_seurat)

oligos_seurat_merged <- Merge_Seurat_List(oligos_seurat_list)

########################################

oligos_seurat_merged <- NormalizeData(oligos_seurat_merged)
oligos_seurat_merged <- FindVariableFeatures(oligos_seurat_merged)

oligos_sketched <- SketchData(object = oligos_seurat_merged, ncells = 30000, method = "LeverageScore", sketched.assay = "sketch")

keep(oligos_sketched)
keep(oligos_sketched, sure = T)

########################################

DefaultAssay(oligos_sketched) <- "sketch" 
oligos_sketched <- FindVariableFeatures(oligos_sketched)
oligos_sketched <- ScaleData(oligos_sketched, vars.to.regress = "pct.mito")
oligos_sketched <- RunPCA(oligos_sketched)
oligos_sketched[["sketch"]] <- JoinLayers(oligos_sketched[["sketch"]])
oligos_sketched <- RunHarmony(object = oligos_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"), 
                             reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(oligos_sketched, ndims = 50, reduction = "harmony")
oligos_sketched <- FindNeighbors(oligos_sketched, reduction = "harmony", dims = 1:15)
oligos_sketched <- RunUMAP(oligos_sketched, reduction = "harmony", dims = 1:15, return.model = T)
oligos_sketched <- FindClusters(oligos_sketched, resolution = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5))

########################################

res_list <- c("0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5")

oligos_markers_list <- list()
oligos_avgexp_list <- list()

for (res in res_list){
  print(paste0("Finding markers and avgexp for res", res))
  
  Idents(oligos_sketched) <- paste0("sketch_snn_res.", res)
  markers <- FindAllMarkers(oligos_sketched, logfc.threshold = 0.5, min.pct = 0.2, only.pos = T)
  avgexp <- as.data.frame(AverageExpression(oligos_sketched, assay = "RNA", layer = "counts"))
  
  oligos_markers_list[[res]] <- markers
  oligos_avgexp_list[[res]] <- avgexp
}

oligos_sketched_meta <- oligos_sketched@meta.data

saveRDS(oligos_markers_list, file = "./analysis/oligodendrocytes/seurat_objects/processing_outs/oligos_sketched_markers_list.rds")
saveRDS(oligos_avgexp_list, file = "./analysis/oligodendrocytes/seurat_objects/processing_outs/oligos_sketched_avgexp_list.rds")
saveRDS(oligos_sketched_meta, file = "./analysis/oligodendrocytes/seurat_objects/processing_outs/oligos_sketched_metadata.rds")

saveRDS(oligos_sketched, file = "./analysis/oligodendrocytes/seurat_objects/oligos_sketched_not_projected.rds")

########################################
########################################
########################################

oligos_sketched <- readRDS("./analysis/oligodendrocytes/seurat_objects/oligos_sketched_not_projected.rds")
oligos_markers_list <- readRDS("./analysis/oligodendrocytes/seurat_objects/processing_outs/oligos_sketched_markers_list.rds")
oligos_sketched_meta <- readRDS("./analysis/oligodendrocytes/seurat_objects/processing_outs/oligos_sketched_metadata.rds")


oligos_sketched_meta <- na.omit(oligos_sketched_meta)
p1 <- clustree(oligos_sketched_meta, prefix = "sketch_snn_res.")


DimPlot_scCustom(oligos_sketched, group.by = "sketch_snn_res.0.1")
DimPlot_scCustom(oligos_sketched, group.by = "sketch_snn_res.0.2")
DimPlot_scCustom(oligos_sketched, group.by = "sketch_snn_res.0.2", split.by = "region")

DotPlot_scCustom(oligos_sketched, group.by = "sketch_snn_res.0.2", 
                 features = c("MOG", "MAG", "PLP1", "OLIG2", "OPALIN"))

########################################

oligos_sketched <- readRDS("./analysis/oligodendrocytes/seurat_objects/oligos_sketched_not_projected.rds")

# projecting oligos back in at res0.2

set.seed(420)

oligos_sketched$dataset_region <- paste0(oligos_sketched$dataset, "_", oligos_sketched$region)

# need to add a metadata col named after RNA layers so that can split "sketch" assay and re-project properly
names(oligos_sketched@assays$RNA@layers)
oligos_sketched@meta.data <- oligos_sketched@meta.data %>%
  mutate(split_sketch = case_when(
    dataset_region == "AMP-PD_DMV" ~ "1.1.1.1.1.1",
    dataset_region == "AMP-PD_GPi" ~ "2.1.1",
    dataset_region == "AMP-PD_M1" ~ "2.1.1.1",
    dataset_region == "AMP-PD_PFC" ~ "2.1.1.1.1",
    dataset_region == "AMP-PD_V1" ~ "2.1.1.1.1.1",
    dataset_region == "Mathys_AnG" ~ "AnG.2.1",
    dataset_region == "Mathys_EC" ~ "EC.2.1",
    dataset_region == "Mathys_HC" ~ "HC.2.1",
    dataset_region == "Mathys_MTG" ~ "MTG.2.1",
    dataset_region == "Mathys_PFC" ~ "PFC.2.1",
    dataset_region == "Mathys_TH" ~ "TH.2.1",
    dataset_region == "Pineda_M1" ~ "M1.2",
    dataset_region == "Pineda_PFC" ~ "PFC.2",
  ))

oligos_sketched[["sketch"]] <- split(oligos_sketched[["sketch"]], f = oligos_sketched$split_sketch)

oligos_sketched <- ProjectIntegration(object = oligos_sketched, sketched.assay = "sketch", 
                                      assay = "RNA", reduction = "harmony")

options(future.globals.maxSize = 8000 * 1024^2)
oligos_sketched <- ProjectData(object = oligos_sketched, sketched.assay = "sketch", assay = "RNA",
                               sketched.reduction = "harmony.full",
                               full.reduction = "harmony.full", dims = 1:15,
                               refdata = list(celltype.full = "sketch_snn_res.0.2"))

oligos_sketched <- RunUMAP(oligos_sketched, reduction = "harmony.full", dims = 1:15,
                           reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(oligos_sketched) <- "RNA"
Idents(oligos_sketched) <- "celltype.full"

saveRDS(oligos_sketched, file = "./analysis/oligodendrocytes/seurat_objects/oligos_sketched_clustered_projected_res_0.2.rds")

########################################
########################################
########################################

oligos_sketched <- readRDS("./analysis/oligodendrocytes/seurat_objects/oligos_sketched_clustered_projected_res_0.2.rds")

oligos_sketched$donor_region <- paste0(oligos_sketched$donor, "_", oligos_sketched$region)

# remove donors w/ <50 cells
t <- oligos_sketched@meta.data %>%
  group_by(donor_region) %>%
  summarise(n = n())

donors_keep <- t$donor_region[t$n > 50]


# remove cluster 6 -- damaged cells
cells_keep <- colnames(oligos_sketched)[!oligos_sketched$celltype.full == "6" & oligos_sketched$donor_region %in% donors_keep]

oligos_filtered <- subset(oligos_sketched, cells = cells_keep)

rm(oligos_sketched)

# re-sketch
oligos_filtered <- FindVariableFeatures(oligos_filtered)

oligos_sketched <- SketchData(object = oligos_filtered, ncells = 30000, method = "LeverageScore", sketched.assay = "sketch")

keep(oligos_sketched)
keep(oligos_sketched, sure = T)

########################################

DefaultAssay(oligos_sketched) <- "sketch" 
oligos_sketched <- FindVariableFeatures(oligos_sketched)
oligos_sketched <- ScaleData(oligos_sketched, vars.to.regress = "pct.mito")
oligos_sketched <- RunPCA(oligos_sketched)
oligos_sketched[["sketch"]] <- JoinLayers(oligos_sketched[["sketch"]])
oligos_sketched <- RunHarmony(object = oligos_sketched, reduction = "pca", group.by.vars = c("dataset", "donor"), 
                              reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "sketch")
ElbowPlot(oligos_sketched, ndims = 50, reduction = "harmony")
oligos_sketched <- FindNeighbors(oligos_sketched, reduction = "harmony", dims = 1:15)
oligos_sketched <- RunUMAP(oligos_sketched, reduction = "harmony", dims = 1:15, return.model = T)
oligos_sketched <- FindClusters(oligos_sketched, resolution = c(0.1, 0.15, 0.2, 0.25))

########################################

DimPlot_scCustom(oligos_sketched, group.by = "sketch_snn_res.0.1")
DimPlot_scCustom(oligos_sketched, group.by = "sketch_snn_res.0.15")
DimPlot_scCustom(oligos_sketched, group.by = "sketch_snn_res.0.2")

DotPlot_scCustom(oligos_sketched, group.by = "sketch_snn_res.0.15", 
                 features = c("MOG", "MAG", "PLP1", "OLIG2", "OPALIN", "RBFOX1", "HSPH1", "CNTN1", "KCNIP4"))

# going w/ res 0.15

########################################

set.seed(420)

oligos_sketched$dataset_region <- paste0(oligos_sketched$dataset, "_", oligos_sketched$region)

# need to add a metadata col named after RNA layers so that can split "sketch" assay and re-project properly
names(oligos_sketched@assays$RNA@layers)
oligos_sketched@meta.data <- oligos_sketched@meta.data %>%
  mutate(split_sketch = case_when(
    dataset_region == "AMP-PD_DMV" ~ "1.1.1.1.1.1",
    dataset_region == "AMP-PD_GPi" ~ "2.1.1",
    dataset_region == "AMP-PD_M1" ~ "2.1.1.1",
    dataset_region == "AMP-PD_PFC" ~ "2.1.1.1.1",
    dataset_region == "AMP-PD_V1" ~ "2.1.1.1.1.1",
    dataset_region == "Mathys_AnG" ~ "AnG.2.1",
    dataset_region == "Mathys_EC" ~ "EC.2.1",
    dataset_region == "Mathys_HC" ~ "HC.2.1",
    dataset_region == "Mathys_MTG" ~ "MTG.2.1",
    dataset_region == "Mathys_PFC" ~ "PFC.2.1",
    dataset_region == "Mathys_TH" ~ "TH.2.1",
    dataset_region == "Pineda_M1" ~ "M1.2",
    dataset_region == "Pineda_PFC" ~ "PFC.2",
  ))

oligos_sketched[["sketch"]] <- split(oligos_sketched[["sketch"]], f = oligos_sketched$split_sketch)

oligos_sketched <- ProjectIntegration(object = oligos_sketched, sketched.assay = "sketch", 
                                      assay = "RNA", reduction = "harmony")

options(future.globals.maxSize = 8000 * 1024^2)
oligos_sketched <- ProjectData(object = oligos_sketched, sketched.assay = "sketch", assay = "RNA",
                               sketched.reduction = "harmony.full",
                               full.reduction = "harmony.full", dims = 1:15,
                               refdata = list(celltype.full = "sketch_snn_res.0.15"))

oligos_sketched <- RunUMAP(oligos_sketched, reduction = "harmony.full", dims = 1:15,
                           reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(oligos_sketched) <- "RNA"
Idents(oligos_sketched) <- "celltype.full"

saveRDS(oligos_sketched, file = "./analysis/oligodendrocytes/seurat_objects/oligos_sketched_filtered_clustered_projected_res_0.15.rds")

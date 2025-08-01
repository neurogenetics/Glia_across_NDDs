library(Seurat)
library(scCustomize)
library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs/")

########################################
########################################
########################################

cortical <- readRDS("./analysis/astrocytes/seurat_objects/cortical_astrocytes_sketched_filtered_clustered_projected_0.15.rds")

########################################

cortical[["RNA"]] <- JoinLayers(cortical[["RNA"]])

#cortical_markers <- FindAllMarkers(cortical, group.by = "celltype.full", assay = "RNA", only.pos = T, logfc.threshold = 0.5, min.pct = 0.4)


# DotPlot_scCustom(cortical, features = c("AQP4", "GFAP", "CABLES1", "SLC38A1", "WDR49", 
#                                         "GRM3",  "CD44", "GRIA1", "CHI3L1", 
#                                        "HSPH1", "DNAJB1",
#                                         "MAP1B", "MYRIP", "EGF", "ABCC1", "APP", "TXNRD1"))

# FeaturePlot_scCustom(cortical, features = "TXNRD1", label = T)
# VlnPlot_scCustom(cortical, features = "GFAP", pt.size = 0)


# 0 - Astro_Protoplasmic_GRM3
# 1 - Astro_Fibrous_CD44
# 2 - Astro_Fibrous_GRIA1
# 3 - Astro_Reactive_CHI3L1
# 4 - Astro_Stress_HSPH1
# 5 - Astro_Stress_TXNRD1 # these express elements of NRF2 signaling pathway (oxidative stress)


Idents(cortical) <- "celltype.full"

cortical <- RenameIdents(cortical, 
                         "0" = "Astro_Protoplasmic_GRM3", 
                         "1" = "Astro_Fibrous_CD44", 
                         "2" = "Astro_Fibrous_GRIA1", 
                         "3" = "Astro_Reactive_SERPINA3", 
                         "4" = "Astro_Stress_HSPH1", 
                         "5" = "Astro_Stress_TXNRD1")

cortical$cluster_anno <- Idents(cortical)

Idents(cortical) <- "cluster_anno"

#DimPlot_scCustom(cortical, label = T, label.box = T)


saveRDS(cortical, file = "./analysis/astrocytes/seurat_objects/cortical_astrocytes_clustered_ANNOTATED.rds")

########################################

# save updated metadata / counts matrices after filtering

cortical_meta <- cortical@meta.data
saveRDS(cortical_meta, file = "./analysis/astrocytes/cortical_astros_metadata_ANNOTATED.rds")

cortical_cts <- LayerData(cortical, layer = "counts", assay = "RNA")
saveRDS(cortical_cts, file = "./analysis/astrocytes/differential_expression/single_cell_cts_matrices/cortical_astros_filtered_sc_cts.rds")

########################################

# get all markers for FGSEA, etc.

markers_all <- FindAllMarkers(cortical, group.by = "cluster_anno", assay = "RNA", only.pos = F, logfc.threshold = 0, min.pct = 0.2, return.thresh = 1)
saveRDS(markers_all, file = "./analysis/astrocytes/cluster_characterization/cortical_astros_markers_ALL_ANNOTATED.rds")

########################################
########################################
########################################

subcortical <- readRDS("./analysis/astrocytes/seurat_objects/subcortical_astrocytes_sketched_filtered_clustered_projected_0.2.rds")

########################################

subcortical[["RNA"]] <- JoinLayers(subcortical[["RNA"]])

# subcortical_markers <- FindAllMarkers(subcortical, group.by = "celltype.full", assay = "RNA", only.pos = T, logfc.threshold = 0.5, min.pct = 0.4)

# t <- subcortical_markers[subcortical_markers$cluster == "8", ]


# DotPlot_scCustom(subcortical, features = c("AQP4", "GFAP", "CABLES1", "SLC38A1", "CD44",
#                                            "SLC6A1", "FRMD4A", "ARHGAP24", "CACNA2D3",
#                                            "L3MBTL4", "ROBO2", "LINC00609", "ADAMTSL3",
#                                            "NRXN3", "CTNNA3", "KSR2",
#                                            "MDGA2",
#                                            "MECOM", "GRIN2A", "NLGN1", "DNAH11",
#                                           "SERPINA3", "CHI3L1", "GRIA2", "EPHA6", "HSPH1", "SLIT2")) + RotatedAxis() 


# VlnPlot_scCustom(subcortical, features = c("CD44", "GFAP"), pt.size = 0)


Idents(subcortical) <- "celltype.full"

subcortical <- RenameIdents(subcortical, 
                         "0" = "Astro_Protoplasmic_CABLES1", 
                         "1" = "Astro_Fibrous_L4MBTL4", 
                         "2" = "Astro_Fibrous_NRXN3", 
                         "3" = "Astro_MDGA2", 
                         "4" = "Astro_Epen_DNAH11", 
                         "5" = "Astro_Reactive_SERPINA3",
                         "6" = "Astro_Protoplasmic_EPHA6",
                         "7" = "Astro_Stress_HSPH1",
                         "8" = "Astro_SLIT2")

subcortical$cluster_anno <- Idents(subcortical)

Idents(subcortical) <- "cluster_anno"

# DimPlot_scCustom(subcortical, label = T, label.box = T)

saveRDS(subcortical, file = "./analysis/astrocytes/seurat_objects/subcortical_astrocytes_clustered_ANNOTATED.rds")

########################################

# save updated metadata / counts matrices after filtering

subcortical_meta <- subcortical@meta.data
saveRDS(subcortical_meta, file = "./analysis/astrocytes/subcortical_astros_metadata_ANNOTATED.rds")

subcortical_cts <- LayerData(subcortical, layer = "counts", assay = "RNA")
saveRDS(subcortical_cts, file = "./analysis/astrocytes/differential_expression/single_cell_cts_matrices/subcortical_astros_filtered_sc_cts.rds")

########################################

# get all markers for FGSEA, etc.

markers_all <- FindAllMarkers(subcortical, group.by = "cluster_anno", assay = "RNA", only.pos = F, logfc.threshold = 0, min.pct = 0.2, return.thresh = 1)
saveRDS(markers_all, file = "./analysis/astrocytes/cluster_characterization/subcortical_astros_markers_ALL_ANNOTATED.rds")

########################################
########################################
########################################

astros_all <- readRDS("./analysis/astrocytes/seurat_objects/all_astros_merged_post_annotation_clustered_projected.rds")

########################################

astros_all[["RNA"]] <- JoinLayers(astros_all[["RNA"]])

markers <- FindAllMarkers(astros_all, group.by = "celltype.full", assay = "RNA", only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)

Idents(astros_all) <- "celltype.full"

DotPlot_scCustom(astros_all, features = c("AQP4", "GFAP", "CABLES1", "SLC38A1", "CD44", "KCND2", "LAMA2")) + RotatedAxis()

astros_all <- RenameIdents(astros_all, 
                            "0" = "Astro_Protoplasmic_GRM3", 
                            "1" = "Astro_Fibrous_DLCK1", 
                            "2" = "Astro_Fibrous_GRIA1", 
                            "3" = "Astro_KCND2", 
                            "4" = "Astro_Reactive_SERPINA3", 
                            "5" = "Astro_LAMA2",
                            "6" = "Astro_Stress_HSPH1")

astros_all$cluster_anno <- Idents(astros_all)

Idents(astros_all) <- "cluster_anno"

DimPlot_scCustom(astros_all, label = T, label.box = T)

saveRDS(astros_all, file = "./analysis/astrocytes/seurat_objects/all_astrocytes_clustered_ANNOTATED.rds")

########################################

# save updated metadata / counts matrices after filtering

all_astros_meta <- astros_all@meta.data
saveRDS(all_astros_meta, file = "./analysis/astrocytes/all_astros_metadata_ANNOTATED.rds")

########################################

# get all markers for FGSEA, etc.

markers_all <- FindAllMarkers(astros_all, group.by = "cluster_anno", assay = "RNA", only.pos = F, logfc.threshold = 0, min.pct = 0.2, return.thresh = 1)
saveRDS(markers_all, file = "./analysis/astrocytes/cluster_characterization/all_astros_markers_ALL_ANNOTATED.rds")

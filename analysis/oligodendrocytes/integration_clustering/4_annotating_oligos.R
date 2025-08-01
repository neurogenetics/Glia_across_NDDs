library(Seurat)
library(scCustomize)
library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs/")

########################################

oligos <- readRDS("./analysis/oligodendrocytes/seurat_objects/oligos_sketched_filtered_clustered_projected_res_0.15.rds")

Idents(oligos) <- "celltype.full"

oligos[["sketch"]] <- JoinLayers(oligos[["sketch"]])
# oligo_markers <- FindAllMarkers(oligos, group.by = "sketch_snn_res.0.15", assay = "sketch", only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
# 
# VlnPlot_scCustom(oligos, features = c("OPALIN"), pt.size = 0)
# DotPlot_scCustom(oligos, features = unique(c("OPALIN", "RBFOX1", "RASGRF1", "ACTN2", "CNTN1", "FRY", "KCNIP4",
#                                              "CSMD1", "SPARCL1", "SYT1", "FGFR1", "BCL6", "NAV2", "HSPH1")))


Idents(oligos) <- "celltype.full"

oligos <- RenameIdents(oligos, 
                         "0" = "Oligo_OPALIN", 
                         "1" = "Oligo_RBFOX1", 
                         "2" = "Oligo_Stress_HSPH1", 
                         "3" = "Oligo_CNTN1", 
                         "4" = "Oligo_CSMD1")

oligos$cluster_anno <- Idents(oligos)

Idents(oligos) <- "cluster_anno"

saveRDS(oligos, file = "./analysis/oligodendrocytes/seurat_objects/oligos_clustered_ANNOTATED.rds")

########################################

# save updated metadata / counts matrices after filtering

oligos_meta <- oligos@meta.data
saveRDS(oligos_meta, file = "./analysis/oligodendrocytes/oligos_celllevel_metadata_ANNOTATED.rds")

########################################

# get all markers for FGSEA, etc.

markers_all <- FindAllMarkers(oligos, group.by = "cluster_anno", assay = "sketch", only.pos = F, logfc.threshold = 0, min.pct = 0.2, return.thresh = 1)
saveRDS(markers_all, file = "./analysis/oligodendrocytes/cluster_characterization/oligos_markers_ALL_ANNOTATED.rds")

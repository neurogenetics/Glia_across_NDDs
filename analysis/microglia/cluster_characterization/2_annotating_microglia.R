library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

########################################

markers <- readRDS("./combined_data/testing_cts_filtering/markers_avgexp/micro_markers_across_resolutions_NO_ctsfilter.rds")
metadata <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata.rds")

########################################

markers_t <- markers[["0.25"]][markers[["0.25"]]$cluster == "12",]

# 0 - Micro_Homeo
# 1 - Micro_DAM_Int1
# 2 - Micro_Inflamm_Stress (== DIMs)
# 3 - BAM
# 4 - Micro_DAM_GPNMB
# 5 - Micro_Inflamm_PCDH9
# 6 - Micro_Phago_CD163
# 7 - Micro_DAM_SPP1
# 8 - Micro_Prolif
# 9 - Micro_DAM_Int2
# 10 - Micro_Inflamm_CD83
# 11 - Monocyte
# 12 - Micro_Inflamm_PCDH9
# 13 - Micro_IFN
# 14 - Lymphocyte

metadata <- metadata %>%
  mutate(cluster_anno = case_when(
    celltype.full == "0" ~ "Micro_Homeo",
    celltype.full == "1" ~ "Micro_DAM_Int1",
    celltype.full == "2" ~ "Micro_Inflamm_Stress",
    celltype.full == "3" ~ "BAM",
    celltype.full == "4" ~ "Micro_DAM_GPNMB",
    celltype.full == "5" ~ "Micro_Inflamm_PCDH9",
    celltype.full == "6" ~ "Micro_Phago_CD163",
    celltype.full == "7" ~ "Micro_DAM_SPP1",
    celltype.full == "8" ~ "Micro_Prolif",
    celltype.full == "9" ~ "Micro_DAM_Int2",
    celltype.full == "10" ~ "Micro_Inflamm_CD83",
    celltype.full == "11" ~ "Monocyte",
    celltype.full == "12" ~ "Micro_Inflamm_PCDH9",
    celltype.full == "13" ~ "Micro_IFN",
    celltype.full == "14" ~ "Lymphocyte",
  ))

saveRDS(metadata, file = "./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")

########################################
########################################
########################################

# need to re-run markers w/ annotations

micro <- readRDS("./combined_data/final_objects/micro_sketched_clustered_projected_NO_ctsfilter_res_0.25.rds")

micro@meta.data <- micro@meta.data %>%
  mutate(cluster_anno = case_when(
    celltype.full == "0" ~ "Micro_Homeo",
    celltype.full == "1" ~ "Micro_DAM_Int1",
    celltype.full == "2" ~ "Micro_Inflamm_Stress",
    celltype.full == "3" ~ "BAM",
    celltype.full == "4" ~ "Micro_DAM_GPNMB",
    celltype.full == "5" ~ "Micro_Inflamm_PCDH9",
    celltype.full == "6" ~ "Micro_Phago_CD163",
    celltype.full == "7" ~ "Micro_DAM_SPP1",
    celltype.full == "8" ~ "Micro_Prolif",
    celltype.full == "9" ~ "Micro_DAM_Int2",
    celltype.full == "10" ~ "Micro_Inflamm_CD83",
    celltype.full == "11" ~ "Monocyte",
    celltype.full == "12" ~ "Micro_Inflamm_PCDH9",
    celltype.full == "13" ~ "Micro_IFN",
    celltype.full == "14" ~ "Lymphocyte",
  ))

Idents(micro) <- "cluster_anno"

micro[["RNA"]] <- JoinLayers(micro[["RNA"]])

markers_anno <- FindAllMarkers(micro, group.by = "cluster_anno", assay = "RNA", only.pos = F, 
                               logfc.threshold = 0, min.pct = 0.2, return.thresh = 1)


saveRDS(markers_anno, file = "./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea_ANNOTATED.rds")
write.csv(markers_anno, file = "./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea_ANNOTATED.csv")

saveRDS(micro, file = "./analysis/microglia/seurat_objects/microglia_annotated.rds")

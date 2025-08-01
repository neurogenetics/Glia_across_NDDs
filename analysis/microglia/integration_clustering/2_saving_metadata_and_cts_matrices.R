library(tidyverse)
library(Seurat)
library(scCustomize)
library(BPCells)

setwd("/data/ADRD/glia_across_NDDs")

########################################

micro <- readRDS("./combined_data/final_objects/micro_sketched_clustered_projected_NO_ctsfilter_res_0.25.rds")

micro_celllevel_meta <- micro@meta.data %>%
  dplyr::select(-orig.ident) %>%
  dplyr::select(-dataset_region)

# edit region names in metadata
micro_celllevel_meta <- micro_celllevel_meta %>%
  mutate(region = case_when(
    region == "PFC" ~ "FC/PFC",
    region == "FC" ~ "FC/PFC",
    region == "TC" ~ "TC/MTG",
    region == "MTG" ~ "TC/MTG",
    region == "OC" ~ "OC/V1", 
    region == "V1" ~ "OC/V1", 
    region == "HC" ~ "HIP",
    TRUE ~ region))

saveRDS(micro_celllevel_meta, file = "./analysis/microglia/all_microglia_celllevel_metadata.rds")

########################################

# extract and save single-cell counts matrices from the different datasets

PD_microglia <- readRDS("./combined_data/processed_celltypes_by_dataset/PD_microglia_allregions_filtered.rds")
AD_microglia <- readRDS("./combined_data/processed_celltypes_by_dataset/AD_microglia_allregions_filtered.rds")
ALS_microglia <- readRDS("./combined_data/processed_celltypes_by_dataset/ALS_microglia_allregions_filtered.rds")
FTD_microglia <- readRDS("./combined_data/processed_celltypes_by_dataset/FTD_microglia_allregions_filtered.rds")

PD_counts <- GetAssayData(PD_microglia, assay = "RNA", layer = "counts")
AD_counts <- GetAssayData(AD_microglia, assay = "RNA", layer = "counts")
ALS_counts <- GetAssayData(ALS_microglia, assay = "RNA", layer = "counts")
FTD_counts <- GetAssayData(FTD_microglia, assay = "RNA", layer = "counts")

saveRDS(PD_counts, file = "./analysis/microglia/differential_expression/single_cell_cts_matrices/AMP-PD_microglia_sc_cts.rds")
saveRDS(AD_counts, file = "./analysis/microglia/differential_expression/single_cell_cts_matrices/Mathys_microglia_sc_cts.rds")
saveRDS(ALS_counts, file = "./analysis/microglia/differential_expression/single_cell_cts_matrices/Pineda_microglia_sc_cts.rds")
saveRDS(FTD_counts, file = "./analysis/microglia/differential_expression/single_cell_cts_matrices/Gerrits_microglia_sc_cts.rds")

########################################

# pseudobulk counts matrices donor x region

AMPPD_cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/AMP-PD_microglia_sc_cts.rds")
Gerrits_cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Gerrits_microglia_sc_cts.rds")
Pineda_cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Pineda_microglia_sc_cts.rds")
Mathys_cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Mathys_microglia_sc_cts.rds")

cts_merged <- cbind(AMPPD_cts, Gerrits_cts, Pineda_cts, Mathys_cts)

celllevel_meta <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")
celllevel_meta$donor_region <- paste0(celllevel_meta$donor, "_", celllevel_meta$region)


pb_cts <- aggregate(as.matrix(t(cts_merged)), by = list(celllevel_meta$donor_region), FUN = sum)
colnames(pb_cts)[1] <- "donor_region"
pb_cts <- pb_cts %>%
  column_to_rownames(var = "donor_region") %>%
  t(.) %>%
  as.data.frame(.)

saveRDS(pb_cts, file = "./analysis/microglia/differential_expression/cts_pseudobulked_donor_x_region.rds")

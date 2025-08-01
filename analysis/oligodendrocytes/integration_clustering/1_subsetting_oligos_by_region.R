library(Seurat)
library(tidyverse)
library(scCustomize)

setwd("/data/ADRD/glia_across_NDDs")

########################################

PD_DMV_merged <- readRDS("./combined_data/processed_regions/pd/DMV_all_preprocessed.rds")
PD_GPi_merged <- readRDS("./combined_data/processed_regions/pd/GPi_all_preprocessed.rds")
PD_M1_merged <- readRDS("./combined_data/processed_regions/pd/M1_all_preprocessed.rds")
PD_PFC_merged <- readRDS("./combined_data/processed_regions/pd/PFC_all_preprocessed.rds")
PD_V1_merged <- readRDS("./combined_data/processed_regions/pd/V1_all_preprocessed.rds")


DimPlot_scCustom(PD_V1_merged)
DotPlot_scCustom(PD_V1_merged, features = c("OPALIN", "MOG", "ST18"))
FeaturePlot_scCustom(PD_V1_merged, features = c("OPALIN", "MOG", "ST18"))
table(PD_V1_merged$harmony_clusters)


PD_DMV_oligos <- subset(PD_DMV_merged, idents = c("0", "1"))
PD_GPi_oligos <- subset(PD_GPi_merged, idents = c("0", "1"))
PD_M1_oligos <- subset(PD_M1_merged, idents = c("0", "1"))
PD_PFC_oligos <- subset(PD_PFC_merged, idents = c("0", "1"))
PD_V1_oligos <- subset(PD_V1_merged, idents = c("0"))


saveRDS(PD_DMV_oligos, file = "./oligos_data/subset_by_dataset_region/PD_DMV_oligos.rds")
saveRDS(PD_GPi_oligos, file = "./oligos_data/subset_by_dataset_region/PD_GPi_oligos.rds")
saveRDS(PD_M1_oligos, file = "./oligos_data/subset_by_dataset_region/PD_M1_oligos.rds")
saveRDS(PD_PFC_oligos, file = "./oligos_data/subset_by_dataset_region/PD_PFC_oligos.rds")
saveRDS(PD_V1_oligos, file = "./oligos_data/subset_by_dataset_region/PD_V1_oligos.rds")

########################################

AD_TH_merged <- readRDS("./combined_data/processed_regions/ad/TH_all_preprocessed.rds")
AD_HIP_merged <- readRDS("./combined_data/processed_regions/ad/HIP_all_preprocessed.rds")
AD_MTG_merged <- readRDS("./combined_data/processed_regions/ad/MTG_all_preprocessed.rds")
AD_PFC_merged <- readRDS("./combined_data/processed_regions/ad/PFC_all_preprocessed.rds")
AD_EC_merged <- readRDS("./combined_data/processed_regions/ad/EC_all_preprocessed.rds")
AD_AnG_merged <- readRDS("./combined_data/processed_regions/ad/AnG_all_preprocessed.rds")


DotPlot_scCustom(AD_AnG_merged, features = c("OPALIN", "MOG", "ST18"))


AD_TH_oligos <- subset(AD_TH_merged, idents = c("0", "1"))
AD_HIP_oligos <- subset(AD_HIP_merged, idents = c("0"))
AD_MTG_oligos <- subset(AD_MTG_merged, idents = c("0"))
AD_PFC_oligos <- subset(AD_PFC_merged, idents = c("0"))
AD_EC_oligos <- subset(AD_EC_merged, idents = c("0"))
AD_AnG_oligos <- subset(AD_AnG_merged, idents = c("0"))


saveRDS(AD_TH_oligos, file = "./oligos_data/subset_by_dataset_region/AD_TH_oligos.rds")
saveRDS(AD_HIP_oligos, file = "./oligos_data/subset_by_dataset_region/AD_HIP_oligos.rds")
saveRDS(AD_MTG_oligos, file = "./oligos_data/subset_by_dataset_region/AD_MTG_oligos.rds")
saveRDS(AD_PFC_oligos, file = "./oligos_data/subset_by_dataset_region/AD_PFC_oligos.rds")
saveRDS(AD_EC_oligos, file = "./oligos_data/subset_by_dataset_region/AD_EC_oligos.rds")
saveRDS(AD_AnG_oligos, file = "./oligos_data/subset_by_dataset_region/AD_AnG_oligos.rds")

########################################

ALS_M1_merged <- readRDS("./combined_data/processed_regions/als/M1_all_preprocessed.rds")
ALS_PFC_merged <- readRDS("./combined_data/processed_regions/als/PFC_all_preprocessed.rds")


DimPlot_scCustom(ALS_EC_merged)
DotPlot_scCustom(ALS_EC_merged, features = c("OPALIN", "MOG", "ST18"))
FeaturePlot_scCustom(ALS_EC_merged, features = c("OPALIN", "MOG", "ST18"))
table(ALS_EC_merged$harmony_clusters)


ALS_M1_oligos <- subset(ALS_M1_merged, idents = c("0", "1"))
ALS_PFC_oligos <- subset(ALS_PFC_merged, idents = c("0"))


saveRDS(ALS_M1_oligos, file = "./oligos_data/subset_by_dataset_region/ALS_M1_oligos.rds")
saveRDS(ALS_PFC_oligos, file = "./oligos_data/subset_by_dataset_region/ALS_PFC_oligos.rds")


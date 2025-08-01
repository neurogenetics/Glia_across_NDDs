library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

Mathys_2023 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Mathys_2023/PCA_GLM_summary.rds")
Green_2024 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_GLM_summary.rds")
SEAAD_PFC <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_GLM_summary.rds")
Macnair_2025 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_GLM_summary.rds")
Rexach_2024 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/PCA_GLM_summary.rds")
Kamath_2022 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/PCA_GLM_summary.rds")
Petrescu_2025 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_GLM_summary.rds")
Prudencio_2015 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_GLM_summary.rds")
Martirosyan_2024 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/PCA_GLM_summary.rds")
Zelic_2025 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/PCA_GLM_summary.rds")
Wang_2024 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/PCA_GLM_summary.rds")

##################################################

Macnair_2025 <- Macnair_2025 %>%
  dplyr::select(-df)

# FDR correction

rep <- rbind(Mathys_2023, Green_2024, SEAAD_PFC, Macnair_2025, Rexach_2024, Kamath_2022, 
             Petrescu_2025, Prudencio_2015, Martirosyan_2024, Zelic_2025, Wang_2024)

rep$padj <- p.adjust(rep$`Pr(>|t|)`, method = "BH")

rep$Estimate <- abs(rep$Estimate)

rep$`-log10FDR` <- -log10(rep$padj)

rep$test <- c("Green_2024_PFC_CERAD", "Green_2024_PFC_Braak", "Green_2024_PFC_MCI", "Green_2024_PFC_AD_dementia",
              "Mathys_2023_PFC_CERAD", "Mathys_2023_PFC_Braak", "Mathys_2023_PFC_MCI", "Mathys_2023_PFC_AD_dementia",
              "Gabitto_2024_PFC_CERAD", "Gabitto_2024_PFC_Braak", "Gabitto_2024_PFC_AD_dementia",
              "Macnair_2025_CTX_WM_NAWM", "Macnair_2025_CTX_WM_AL", "Macnair_2025_CTX_WM_CAL",
              "Macnair_2025_CTX_WM_CIL", "Macnair_2025_CTX_WM_RL",
              "Macnair_2025_CTX_GM_NAGM", "Macnair_2025_CTX_GM_GML",
              "Rexach_2024_INS_AD", "Rexach_2024_INS_PSP", "Rexach_2024_INS_PiD",
              "Kamath_2022_SN_PD/DLB",
              "Petrescu_2025_BA44_ALS", "Petrescu_2025_BA44_ALS_borderline", "Petrescu_2025_BA44_ALSci",
              "Petrescu_2025_BA46_ALS", "Petrescu_2025_BA46_ALS_borderline", "Petrescu_2025_BA46_ALSci",
              "Prudencio_2015_FC_sALS", "Prudencio_2015_FC_C9ALS",
              "Prudencio_2015_CB_sALS", "Prudencio_2015_CB_C9ALS",
              "Martirosyan_2024_SN_PD",
              "Zelic_2025_SpC_ALS",
              "Wang_2024_SN_PD/PDD")

saveRDS(rep, "./analysis/microglia/DAM_signature_replication/GLM_summary_merged_FDRcorr.rds")

library(tidyverse)
library(DESeq2)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################
##################################################
##################################################

# Green 2024

##################################################

# CERAD

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$individualID]

pb_cts <- pb_cts[, meta$individualID]

meta <- meta %>%
  column_to_rownames(var = "individualID")

meta$pmi_scaled <- scale(meta$pmi)

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ cerad_numerical + age_binned + pmi_scaled + msex + batch)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Green_2024_cerad_deseq.rds")

##################################################

# Braak

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$individualID]

pb_cts <- pb_cts[, meta$individualID]

meta <- meta %>%
  column_to_rownames(var = "individualID")

meta$pmi_scaled <- scale(meta$pmi)

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ braak_scaled + age_binned + pmi_scaled + msex + batch)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Green_2024_braak_deseq.rds")

##################################################

# cognitive status

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$individualID]

pb_cts <- pb_cts[, meta$individualID]

meta <- meta %>%
  column_to_rownames(var = "individualID")

meta$cogdx <- factor(meta$cogdx, levels = c("HC", "MCI", "AD_dementia"))

meta$pmi_scaled <- scale(meta$pmi)

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ cogdx + age_binned + pmi_scaled + msex + batch)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Green_2024_cogdx_deseq.rds")

##################################################
##################################################
##################################################

# Mathys 2023

##################################################

# CERAD

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Mathys_2023/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Mathys_2023/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$projid]

pb_cts <- pb_cts[, meta$projid]

meta <- meta %>%
  column_to_rownames(var = "projid")

meta$pmi_scaled <- scale(meta$pmi)

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ cerad_numerical + age_binned + pmi_scaled + msex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Mathys_2023_cerad_deseq.rds")

##################################################

# Braak

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Mathys_2023/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Mathys_2023/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$projid]

pb_cts <- pb_cts[, meta$projid]

meta <- meta %>%
  column_to_rownames(var = "projid")

meta$pmi_scaled <- scale(meta$pmi)

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ braak_scaled + age_binned + pmi_scaled + msex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Mathys_2023_braak_deseq.rds")

##################################################

# cognitive status

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Mathys_2023/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Mathys_2023/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$projid]

pb_cts <- pb_cts[, meta$projid]

meta <- meta %>%
  column_to_rownames(var = "projid")

meta$cogdx <- factor(meta$cogdx, levels = c("HC", "MCI", "AD_dementia"))

meta$pmi_scaled <- scale(meta$pmi)

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ cogdx + age_binned + pmi_scaled + msex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Mathys_2023_cogdx_deseq.rds")

##################################################
##################################################
##################################################

# SEAAD PFC

##################################################

# CERAD

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$donor]

pb_cts <- pb_cts[, meta$donor]

pb_cts <- round(pb_cts)

meta <- meta %>%
  column_to_rownames(var = "donor")

meta$pmi_scaled <- scale(meta$pmi)
meta$age_scaled <- scale(meta$age)

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ cerad_numerical + age_scaled + pmi_scaled + sex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/SEAAD_PFC_cerad_deseq.rds")

##################################################

# Braak

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$donor]

pb_cts <- pb_cts[, meta$donor]

pb_cts <- round(pb_cts)

meta <- meta %>%
  column_to_rownames(var = "donor")

meta$pmi_scaled <- scale(meta$pmi)
meta$age_scaled <- scale(meta$age)

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ braak_scaled + age_scaled + pmi_scaled + sex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/SEAAD_PFC_braak_deseq.rds")

##################################################

# cognitive status

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$donor]

pb_cts <- pb_cts[, meta$donor]

pb_cts <- round(pb_cts)

meta <- meta %>%
  column_to_rownames(var = "donor")

meta$pmi_scaled <- scale(meta$pmi)
meta$age_scaled <- scale(meta$age)

meta$cogdx <- factor(meta$cogdx, levels = c("HC", "AD_dementia"))

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ cogdx + age_scaled + pmi_scaled + sex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/SEAAD_PFC_cogdx_deseq.rds")

##################################################
##################################################
##################################################

# Rexach 2024

##################################################

# neuropathological diagnosis

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/PCA_DAM_final.rds")

meta$library_id <- gsub("_", ".", meta$library_id)

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$library_id]

pb_cts <- pb_cts[, meta$library_id]

meta <- meta %>%
  column_to_rownames(var = "library_id")

meta$pmi_scaled <- scale(meta$pmi)
meta$age_scaled <- scale(meta$age)

meta$npdx1 <- factor(meta$npdx1, levels = c("HC", "AD", "PSP", "PiD"))

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ npdx1 + age_scaled + pmi_scaled + sex + finalsite)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Rexach_2024_npdx_deseq.rds")

##################################################
##################################################
##################################################

# Macnair 2025

##################################################

# white matter

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_DAM_final_WM.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$sample_id_anon]

pb_cts <- pb_cts[, meta$sample_id_anon]

meta <- meta %>%
  column_to_rownames(var = "sample_id_anon")

meta$pmi_scaled <- scale(meta$pmi)
meta$age_scaled <- scale(meta$age)

meta$lesion_type <- factor(meta$lesion_type, levels = c("WM", "NAWM", "AL", "CAL", "CIL", "RL"))

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ lesion_type + age_scaled + pmi_scaled + sex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Macnair_2025_WM_deseq.rds")

##################################################

# grey matter

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_DAM_final_GM.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$sample_id_anon]

pb_cts <- pb_cts[, meta$sample_id_anon]

meta <- meta %>%
  column_to_rownames(var = "sample_id_anon")

meta$pmi_scaled <- scale(meta$pmi)
meta$age_scaled <- scale(meta$age)

meta$lesion_type <- factor(meta$lesion_type, levels = c("GM", "NAGM", "GML"))

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ lesion_type + age_scaled + pmi_scaled + sex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Macnair_2025_GM_deseq.rds")

##################################################
##################################################
##################################################

# Petrescu 2025

##################################################

# BA44

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_DAM_final_BA44.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$library_name]

pb_cts <- pb_cts[, meta$library_name]

meta <- meta %>%
  column_to_rownames(var = "library_name")

meta$age_scaled <- scale(meta$Age)

meta$Classifcation <- factor(meta$Classifcation, levels = c("Control", "ALS", "ALS - borderline", "ALSci"))

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ Classifcation + age_scaled + Sex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Petrescu_2025_BA44_deseq.rds")

##################################################

# BA46

#########################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_DAM_final_BA46.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$library_name]

pb_cts <- pb_cts[, meta$library_name]

meta <- meta %>%
  column_to_rownames(var = "library_name")

meta$age_scaled <- scale(meta$Age)

meta$Classifcation <- factor(meta$Classifcation, levels = c("Control", "ALS", "ALS - borderline", "ALSci"))

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ Classifcation + age_scaled + Sex)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Petrescu_2025_BA46_deseq.rds")

##################################################
##################################################
##################################################

# Wang 2024

##################################################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/PCA_DAM_final.rds")

colnames(pb_cts) <- gsub("_.*", "", colnames(pb_cts))

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$donor]

pb_cts <- pb_cts[, meta$donor]

meta <- meta %>%
  column_to_rownames(var = "donor")

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ group + sex + age_scaled + pmi_scaled)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Wang_2024_dx_PD_deseq.rds")

##################################################
##################################################
##################################################

# Prudencio 2015

##################################################

# frontal cortex

#########################

cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/cts_fcx.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_DAM_final_fcx.rds")

cts <- cts[, colnames(cts) %in% meta$donor]

cts <- cts[, meta$donor]

meta <- meta %>%
  column_to_rownames(var = "donor")

dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = meta, 
                              design = ~ group + sex + age_scaled + pmi_scaled)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Prudencio_2015_fcx_deseq.rds")

##################################################

# cerebellum

#########################

cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/cts_cereb.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_DAM_final_cereb.rds")

cts <- cts[, colnames(cts) %in% meta$donor]

cts <- cts[, meta$donor]

meta <- meta %>%
  column_to_rownames(var = "donor")

dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = meta, 
                              design = ~ group + sex + age_scaled + pmi_scaled)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Prudencio_2015_cereb_deseq.rds")

##################################################
##################################################
##################################################

# Martirosyan 2024

##################################################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$donor]

pb_cts <- pb_cts[, meta$donor]

meta <- meta %>%
  column_to_rownames(var = "donor")

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ group + sex + age_scaled + pmi_scaled)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Martirosyan_2024_dx_PD_deseq.rds")

##################################################
##################################################
##################################################

# Zelic 2025

##################################################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/PCA_DAM_final.rds")

colnames(pb_cts) <- gsub("^[^_]*_", "", colnames(pb_cts))

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$donor]

pb_cts <- pb_cts[, meta$donor]

meta <- meta %>%
  column_to_rownames(var = "donor")

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ group + sex + age_scaled + pmi_scaled)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Zelic_2025_dx_ALS_deseq.rds")

##################################################
##################################################
##################################################

# Kamath 2022

##################################################

pb_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/microglia_pb_cts.rds")
meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/PCA_DAM_final.rds")

pb_cts <- pb_cts[, colnames(pb_cts) %in% meta$donor]

pb_cts <- pb_cts[, meta$donor]

meta <- meta %>%
  column_to_rownames(var = "donor")

dds <- DESeqDataSetFromMatrix(countData = pb_cts, 
                              colData = meta, 
                              design = ~ group + sex + age_scaled + pmi_scaled)

dds <- DESeq(dds)

saveRDS(dds, "./analysis/microglia/DAM_signature_replication/DESeq_results/Kamath_2022_dx_PD_deseq.rds")

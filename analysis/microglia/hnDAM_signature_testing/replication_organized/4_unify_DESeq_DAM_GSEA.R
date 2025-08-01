library(tidyverse)
library(fgsea)
library(DESeq2)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
pathway <- list(DAM_genes)
names(pathway) <- "DAM"

##################################################
##################################################
##################################################

# Rexach 2024

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Rexach_2024_npdx_deseq.rds")

resultsNames(deseq_res)

#########################

# AD

res_AD <- as.data.frame(results(deseq_res, name = "npdx1_AD_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_AD$t <- sign(res_AD$log2FoldChange) * -log10(res_AD$pvalue)
res_AD <- res_AD %>%
  arrange(desc(t))

ranking_AD <- res_AD$t
names(ranking_AD) <- res_AD$gene

gsea_AD <- fgsea(pathways = pathway, stats = ranking_AD, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_AD)

#########################

# PSP

res_PSP <- as.data.frame(results(deseq_res, name = "npdx1_PSP_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_PSP$t <- sign(res_PSP$log2FoldChange) * -log10(res_PSP$pvalue)
res_PSP <- res_PSP %>%
  arrange(desc(t))

ranking_PSP <- res_PSP$t
names(ranking_PSP) <- res_PSP$gene

gsea_PSP <- fgsea(pathways = pathway, stats = ranking_PSP, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_PSP)

#########################

# PiD

res_PiD <- as.data.frame(results(deseq_res, name = "npdx1_PiD_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_PiD$t <- sign(res_PiD$log2FoldChange) * -log10(res_PiD$pvalue)
res_PiD <- res_PiD %>%
  arrange(desc(t))

ranking_PiD <- res_PiD$t
names(ranking_PiD) <- res_PiD$gene

gsea_PiD <- fgsea(pathways = pathway, stats = ranking_PiD, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_PiD)

#########################

# make summary

Rexach_res <- rbind(gsea_AD, gsea_PSP, gsea_PiD)

##################################################
##################################################
##################################################

# Macnair 2025

##################################################

# white matter

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Macnair_2025_WM_deseq.rds")

resultsNames(deseq_res)

#########################

# NAWM

res_NAWM <- as.data.frame(results(deseq_res, name = "lesion_type_NAWM_vs_WM")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_NAWM$t <- sign(res_NAWM$log2FoldChange) * -log10(res_NAWM$pvalue)
res_NAWM <- res_NAWM %>%
  arrange(desc(t))

ranking_NAWM <- res_NAWM$t
names(ranking_NAWM) <- res_NAWM$gene

gsea_NAWM <- fgsea(pathways = pathway, stats = ranking_NAWM, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_NAWM)

#########################

# AL

res_AL <- as.data.frame(results(deseq_res, name = "lesion_type_AL_vs_WM")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_AL$t <- sign(res_AL$log2FoldChange) * -log10(res_AL$pvalue)
res_AL <- res_AL %>%
  arrange(desc(t))

ranking_AL <- res_AL$t
names(ranking_AL) <- res_AL$gene

gsea_AL <- fgsea(pathways = pathway, stats = ranking_AL, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_AL)

#########################

# CAL

res_CAL <- as.data.frame(results(deseq_res, name = "lesion_type_CAL_vs_WM")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_CAL$t <- sign(res_CAL$log2FoldChange) * -log10(res_CAL$pvalue)
res_CAL <- res_CAL %>%
  arrange(desc(t))

ranking_CAL <- res_CAL$t
names(ranking_CAL) <- res_CAL$gene

gsea_CAL <- fgsea(pathways = pathway, stats = ranking_CAL, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_CAL)

#########################

# CIL

res_CIL <- as.data.frame(results(deseq_res, name = "lesion_type_CIL_vs_WM")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_CIL$t <- sign(res_CIL$log2FoldChange) * -log10(res_CIL$pvalue)
res_CIL <- res_CIL %>%
  arrange(desc(t))

ranking_CIL <- res_CIL$t
names(ranking_CIL) <- res_CIL$gene

gsea_CIL <- fgsea(pathways = pathway, stats = ranking_CIL, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_CIL)

#########################

# RL

res_RL <- as.data.frame(results(deseq_res, name = "lesion_type_RL_vs_WM")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_RL$t <- sign(res_RL$log2FoldChange) * -log10(res_RL$pvalue)
res_RL <- res_RL %>%
  arrange(desc(t))

ranking_RL <- res_RL$t
names(ranking_RL) <- res_RL$gene

gsea_RL <- fgsea(pathways = pathway, stats = ranking_RL, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_RL)

##################################################

# grey matter

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Macnair_2025_GM_deseq.rds")

resultsNames(deseq_res)

#########################

# NAGM

res_NAGM <- as.data.frame(results(deseq_res, name = "lesion_type_NAGM_vs_GM")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_NAGM$t <- sign(res_NAGM$log2FoldChange) * -log10(res_NAGM$pvalue)
res_NAGM <- res_NAGM %>%
  arrange(desc(t))

ranking_NAGM <- res_NAGM$t
names(ranking_NAGM) <- res_NAGM$gene

gsea_NAGM <- fgsea(pathways = pathway, stats = ranking_NAGM, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_NAGM)

#########################

# GML

res_GML <- as.data.frame(results(deseq_res, name = "lesion_type_GML_vs_GM")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_GML$t <- sign(res_GML$log2FoldChange) * -log10(res_GML$pvalue)
res_GML <- res_GML %>%
  arrange(desc(t))

ranking_GML <- res_GML$t
names(ranking_GML) <- res_GML$gene

gsea_GML <- fgsea(pathways = pathway, stats = ranking_GML, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_GML)

#########################

# make summary

Macnair_res <- rbind(gsea_NAWM, gsea_AL, gsea_CAL, gsea_CIL, gsea_RL, gsea_NAGM, gsea_GML)

##################################################
##################################################
##################################################

# Green 2024

##################################################

# CERAD

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Green_2024_cerad_deseq.rds")

resultsNames(deseq_res)

res_Green_cerad <- as.data.frame(results(deseq_res, name = "cerad_numerical")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Green_cerad$t <- sign(res_Green_cerad$log2FoldChange) * -log10(res_Green_cerad$pvalue)
res_Green_cerad <- res_Green_cerad %>%
  arrange(desc(t))

ranking_Green_cerad <- res_Green_cerad$t
names(ranking_Green_cerad) <- res_Green_cerad$gene

gsea_Green_cerad <- fgsea(pathways = pathway, stats = ranking_Green_cerad, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_Green_cerad)

#########################

# Braak

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Green_2024_braak_deseq.rds")

resultsNames(deseq_res)

res_Green_braak <- as.data.frame(results(deseq_res, name = "braak_scaled")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Green_braak$t <- sign(res_Green_braak$log2FoldChange) * -log10(res_Green_braak$pvalue)
res_Green_braak <- res_Green_braak %>%
  arrange(desc(t))

ranking_Green_braak <- res_Green_braak$t
names(ranking_Green_braak) <- res_Green_braak$gene

gsea_Green_braak <- fgsea(pathways = pathway, stats = ranking_Green_braak, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_Green_braak)

#########################

# cogdx - MCI

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Green_2024_cogdx_deseq.rds")

resultsNames(deseq_res)

res_Green_MCI <- as.data.frame(results(deseq_res, name = "cogdx_MCI_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Green_MCI$t <- sign(res_Green_MCI$log2FoldChange) * -log10(res_Green_MCI$pvalue)
res_Green_MCI <- res_Green_MCI %>%
  arrange(desc(t))

ranking_Green_MCI <- res_Green_MCI$t
names(ranking_Green_MCI) <- res_Green_MCI$gene

gsea_Green_MCI <- fgsea(pathways = pathway, stats = ranking_Green_MCI, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_Green_MCI)

#########################

# cogdx - AD_dementia

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Green_2024_cogdx_deseq.rds")

resultsNames(deseq_res)

res_Green_AD_dementia <- as.data.frame(results(deseq_res, name = "cogdx_AD_dementia_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Green_AD_dementia$t <- sign(res_Green_AD_dementia$log2FoldChange) * -log10(res_Green_AD_dementia$pvalue)
res_Green_AD_dementia <- res_Green_AD_dementia %>%
  arrange(desc(t))

ranking_Green_AD_dementia <- res_Green_AD_dementia$t
names(ranking_Green_AD_dementia) <- res_Green_AD_dementia$gene

gsea_Green_AD_dementia <- fgsea(pathways = pathway, stats = ranking_Green_AD_dementia, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_Green_AD_dementia)

#########################

# make summary

Green_res <- rbind(gsea_Green_cerad, gsea_Green_braak, gsea_Green_MCI, gsea_Green_AD_dementia)

##################################################
##################################################
##################################################

# Mathys 2023

##################################################

# CERAD

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Mathys_2023_cerad_deseq.rds")

resultsNames(deseq_res)

res_Mathys_cerad <- as.data.frame(results(deseq_res, name = "cerad_numerical")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Mathys_cerad$t <- sign(res_Mathys_cerad$log2FoldChange) * -log10(res_Mathys_cerad$pvalue)
res_Mathys_cerad <- res_Mathys_cerad %>%
  arrange(desc(t))

ranking_Mathys_cerad <- res_Mathys_cerad$t
names(ranking_Mathys_cerad) <- res_Mathys_cerad$gene

gsea_Mathys_cerad <- fgsea(pathways = pathway, stats = ranking_Mathys_cerad, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_Mathys_cerad)

#########################

# Braak

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Mathys_2023_braak_deseq.rds")

resultsNames(deseq_res)

res_Mathys_braak <- as.data.frame(results(deseq_res, name = "braak_scaled")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Mathys_braak$t <- sign(res_Mathys_braak$log2FoldChange) * -log10(res_Mathys_braak$pvalue)
res_Mathys_braak <- res_Mathys_braak %>%
  arrange(desc(t))

ranking_Mathys_braak <- res_Mathys_braak$t
names(ranking_Mathys_braak) <- res_Mathys_braak$gene

gsea_Mathys_braak <- fgsea(pathways = pathway, stats = ranking_Mathys_braak, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_Mathys_braak)

#########################

# cogdx - MCI

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Mathys_2023_cogdx_deseq.rds")

resultsNames(deseq_res)

res_Mathys_MCI <- as.data.frame(results(deseq_res, name = "cogdx_MCI_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Mathys_MCI$t <- sign(res_Mathys_MCI$log2FoldChange) * -log10(res_Mathys_MCI$pvalue)
res_Mathys_MCI <- res_Mathys_MCI %>%
  arrange(desc(t))

ranking_Mathys_MCI <- res_Mathys_MCI$t
names(ranking_Mathys_MCI) <- res_Mathys_MCI$gene

gsea_Mathys_MCI <- fgsea(pathways = pathway, stats = ranking_Mathys_MCI, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_Mathys_MCI)

#########################

# cogdx - AD_dementia

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Mathys_2023_cogdx_deseq.rds")

resultsNames(deseq_res)

res_Mathys_AD_dementia <- as.data.frame(results(deseq_res, name = "cogdx_AD_dementia_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Mathys_AD_dementia$t <- sign(res_Mathys_AD_dementia$log2FoldChange) * -log10(res_Mathys_AD_dementia$pvalue)
res_Mathys_AD_dementia <- res_Mathys_AD_dementia %>%
  arrange(desc(t))

ranking_Mathys_AD_dementia <- res_Mathys_AD_dementia$t
names(ranking_Mathys_AD_dementia) <- res_Mathys_AD_dementia$gene

gsea_Mathys_AD_dementia <- fgsea(pathways = pathway, stats = ranking_Mathys_AD_dementia, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_Mathys_AD_dementia)

#########################

# make summary

Mathys_res <- rbind(gsea_Mathys_cerad, gsea_Mathys_braak, gsea_Mathys_MCI, gsea_Mathys_AD_dementia)

##################################################
##################################################
##################################################

# SEAAD PFC

##################################################

# CERAD

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/SEAAD_PFC_cerad_deseq.rds")

resultsNames(deseq_res)

res_SEAAD_cerad <- as.data.frame(results(deseq_res, name = "cerad_numerical")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_SEAAD_cerad$t <- sign(res_SEAAD_cerad$log2FoldChange) * -log10(res_SEAAD_cerad$pvalue)
res_SEAAD_cerad <- res_SEAAD_cerad %>%
  arrange(desc(t))

ranking_SEAAD_cerad <- res_SEAAD_cerad$t
names(ranking_SEAAD_cerad) <- res_SEAAD_cerad$gene

gsea_SEAAD_cerad <- fgsea(pathways = pathway, stats = ranking_SEAAD_cerad, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_SEAAD_cerad)

#########################

# Braak

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/SEAAD_PFC_braak_deseq.rds")

resultsNames(deseq_res)

res_SEAAD_braak <- as.data.frame(results(deseq_res, name = "braak_scaled")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_SEAAD_braak$t <- sign(res_SEAAD_braak$log2FoldChange) * -log10(res_SEAAD_braak$pvalue)
res_SEAAD_braak <- res_SEAAD_braak %>%
  arrange(desc(t))

ranking_SEAAD_braak <- res_SEAAD_braak$t
names(ranking_SEAAD_braak) <- res_SEAAD_braak$gene

gsea_SEAAD_braak <- fgsea(pathways = pathway, stats = ranking_SEAAD_braak, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_SEAAD_braak)

#########################

# cogdx - AD_dementia

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/SEAAD_PFC_cogdx_deseq.rds")

resultsNames(deseq_res)

res_SEAAD_AD_dementia <- as.data.frame(results(deseq_res, name = "cogdx_AD_dementia_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_SEAAD_AD_dementia$t <- sign(res_SEAAD_AD_dementia$log2FoldChange) * -log10(res_SEAAD_AD_dementia$pvalue)
res_SEAAD_AD_dementia <- res_SEAAD_AD_dementia %>%
  arrange(desc(t))

ranking_SEAAD_AD_dementia <- res_SEAAD_AD_dementia$t
names(ranking_SEAAD_AD_dementia) <- res_SEAAD_AD_dementia$gene

gsea_SEAAD_AD_dementia <- fgsea(pathways = pathway, stats = ranking_SEAAD_AD_dementia, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_SEAAD_AD_dementia)

#########################

# make summary

SEAAD_res <- rbind(gsea_SEAAD_cerad, gsea_SEAAD_braak, gsea_SEAAD_AD_dementia)

##################################################
##################################################
##################################################

# Petrescu 2025

##################################################

# BA44

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Petrescu_2025_BA44_deseq.rds")

resultsNames(deseq_res)

#########################

# ALS

res_ALS_BA44 <- as.data.frame(results(deseq_res, name = "Classifcation_ALS_vs_Control")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_ALS_BA44$t <- sign(res_ALS_BA44$log2FoldChange) * -log10(res_ALS_BA44$pvalue)
res_ALS_BA44 <- res_ALS_BA44 %>%
  arrange(desc(t))

ranking_ALS_BA44 <- res_ALS_BA44$t
names(ranking_ALS_BA44) <- res_ALS_BA44$gene

gsea_ALS_BA44 <- fgsea(pathways = pathway, stats = ranking_ALS_BA44, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_ALS_BA44)

#########################

# ALS - borderline

res_ALS_borderline_BA44 <- as.data.frame(results(deseq_res, name = "Classifcation_ALS...borderline_vs_Control")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_ALS_borderline_BA44$t <- sign(res_ALS_borderline_BA44$log2FoldChange) * -log10(res_ALS_borderline_BA44$pvalue)
res_ALS_borderline_BA44 <- res_ALS_borderline_BA44 %>%
  arrange(desc(t))

ranking_ALS_borderline_BA44 <- res_ALS_borderline_BA44$t
names(ranking_ALS_borderline_BA44) <- res_ALS_borderline_BA44$gene

gsea_ALS_borderline_BA44 <- fgsea(pathways = pathway, stats = ranking_ALS_borderline_BA44, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_ALS_borderline_BA44)

#########################

# ALSci

res_ALSci_BA44 <- as.data.frame(results(deseq_res, name = "Classifcation_ALSci_vs_Control")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_ALSci_BA44$t <- sign(res_ALSci_BA44$log2FoldChange) * -log10(res_ALSci_BA44$pvalue)
res_ALSci_BA44 <- res_ALSci_BA44 %>%
  arrange(desc(t))

ranking_ALSci_BA44 <- res_ALSci_BA44$t
names(ranking_ALSci_BA44) <- res_ALSci_BA44$gene

gsea_ALSci_BA44 <- fgsea(pathways = pathway, stats = ranking_ALSci_BA44, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_ALSci_BA44)

##################################################

# BA46

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Petrescu_2025_BA46_deseq.rds")

resultsNames(deseq_res)

#########################

# ALS

res_ALS_BA46 <- as.data.frame(results(deseq_res, name = "Classifcation_ALS_vs_Control")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_ALS_BA46$t <- sign(res_ALS_BA46$log2FoldChange) * -log10(res_ALS_BA46$pvalue)
res_ALS_BA46 <- res_ALS_BA46 %>%
  arrange(desc(t))

ranking_ALS_BA46 <- res_ALS_BA46$t
names(ranking_ALS_BA46) <- res_ALS_BA46$gene

gsea_ALS_BA46 <- fgsea(pathways = pathway, stats = ranking_ALS_BA46, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_ALS_BA46)

#########################

# ALS - borderline

res_ALS_borderline_BA46 <- as.data.frame(results(deseq_res, name = "Classifcation_ALS...borderline_vs_Control")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_ALS_borderline_BA46$t <- sign(res_ALS_borderline_BA46$log2FoldChange) * -log10(res_ALS_borderline_BA46$pvalue)
res_ALS_borderline_BA46 <- res_ALS_borderline_BA46 %>%
  arrange(desc(t))

ranking_ALS_borderline_BA46 <- res_ALS_borderline_BA46$t
names(ranking_ALS_borderline_BA46) <- res_ALS_borderline_BA46$gene

gsea_ALS_borderline_BA46 <- fgsea(pathways = pathway, stats = ranking_ALS_borderline_BA46, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_ALS_borderline_BA46)

#########################

# ALSci

res_ALSci_BA46 <- as.data.frame(results(deseq_res, name = "Classifcation_ALSci_vs_Control")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_ALSci_BA46$t <- sign(res_ALSci_BA46$log2FoldChange) * -log10(res_ALSci_BA46$pvalue)
res_ALSci_BA46 <- res_ALSci_BA46 %>%
  arrange(desc(t))

ranking_ALSci_BA46 <- res_ALSci_BA46$t
names(ranking_ALSci_BA46) <- res_ALSci_BA46$gene

gsea_ALSci_BA46 <- fgsea(pathways = pathway, stats = ranking_ALSci_BA46, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_ALSci_BA46)

#########################

# make summary

Petrescu_res <- rbind(gsea_ALS_BA44, gsea_ALS_borderline_BA44, gsea_ALSci_BA44,
                      gsea_ALS_BA46, gsea_ALS_borderline_BA46, gsea_ALSci_BA46)

##################################################
##################################################
##################################################

# Wang 2024

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Wang_2024_dx_PD_deseq.rds")

resultsNames(deseq_res)

#########################

# PD

res_PD <- as.data.frame(results(deseq_res, name = "group_PD_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_PD$t <- sign(res_PD$log2FoldChange) * -log10(res_PD$pvalue)
res_PD <- res_PD %>%
  arrange(desc(t))

ranking_PD <- res_PD$t
names(ranking_PD) <- res_PD$gene

gsea_PD <- fgsea(pathways = pathway, stats = ranking_PD, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_PD)

##################################################
##################################################
##################################################

# Prudencio 2015

##################################################

# frontal cortex

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Prudencio_2015_fcx_deseq.rds")

resultsNames(deseq_res)

#########################

# sALS

res_sALS_fcx <- as.data.frame(results(deseq_res, name = "group_sALS_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_sALS_fcx$t <- sign(res_sALS_fcx$log2FoldChange) * -log10(res_sALS_fcx$pvalue)
res_sALS_fcx <- res_sALS_fcx %>%
  arrange(desc(t))

ranking_sALS_fcx <- res_sALS_fcx$t
names(ranking_sALS_fcx) <- res_sALS_fcx$gene

gsea_sALS_fcx <- fgsea(pathways = pathway, stats = ranking_sALS_fcx, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_sALS_fcx)

#########################

# C9ALS

res_C9ALS_fcx <- as.data.frame(results(deseq_res, name = "group_C9ALS_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_C9ALS_fcx$t <- sign(res_C9ALS_fcx$log2FoldChange) * -log10(res_C9ALS_fcx$pvalue)
res_C9ALS_fcx <- res_C9ALS_fcx %>%
  arrange(desc(t))

ranking_C9ALS_fcx <- res_C9ALS_fcx$t
names(ranking_C9ALS_fcx) <- res_C9ALS_fcx$gene

gsea_C9ALS_fcx <- fgsea(pathways = pathway, stats = ranking_C9ALS_fcx, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_C9ALS_fcx)

##################################################

# cerebellum

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Prudencio_2015_cereb_deseq.rds")

resultsNames(deseq_res)

#########################

# sALS

res_sALS_cereb <- as.data.frame(results(deseq_res, name = "group_sALS_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_sALS_cereb$t <- sign(res_sALS_cereb$log2FoldChange) * -log10(res_sALS_cereb$pvalue)
res_sALS_cereb <- res_sALS_cereb %>%
  arrange(desc(t))

ranking_sALS_cereb <- res_sALS_cereb$t
names(ranking_sALS_cereb) <- res_sALS_cereb$gene

gsea_sALS_cereb <- fgsea(pathways = pathway, stats = ranking_sALS_cereb, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_sALS_cereb)

#########################

# C9ALS

res_C9ALS_cereb <- as.data.frame(results(deseq_res, name = "group_C9ALS_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_C9ALS_cereb$t <- sign(res_C9ALS_cereb$log2FoldChange) * -log10(res_C9ALS_cereb$pvalue)
res_C9ALS_cereb <- res_C9ALS_cereb %>%
  arrange(desc(t))

ranking_C9ALS_cereb <- res_C9ALS_cereb$t
names(ranking_C9ALS_cereb) <- res_C9ALS_cereb$gene

gsea_C9ALS_cereb <- fgsea(pathways = pathway, stats = ranking_C9ALS_cereb, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_C9ALS_cereb)

#########################

# make summary

Prudencio_res <- rbind(gsea_sALS_fcx, gsea_C9ALS_fcx,
                       gsea_sALS_cereb, gsea_C9ALS_cereb)

##################################################
##################################################
##################################################

# Martirosyan 2024

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Martirosyan_2024_dx_PD_deseq.rds")

resultsNames(deseq_res)

#########################

# PD

res_PD_Martirosyan <- as.data.frame(results(deseq_res, name = "group_PD_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_PD_Martirosyan$t <- sign(res_PD_Martirosyan$log2FoldChange) * -log10(res_PD_Martirosyan$pvalue)
res_PD_Martirosyan <- res_PD_Martirosyan %>%
  arrange(desc(t))

ranking_PD_Martirosyan <- res_PD_Martirosyan$t
names(ranking_PD_Martirosyan) <- res_PD_Martirosyan$gene

gsea_PD_Martirosyan <- fgsea(pathways = pathway, stats = ranking_PD_Martirosyan, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_PD_Martirosyan)

##################################################
##################################################
##################################################

# Zelic 2025

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Zelic_2025_dx_ALS_deseq.rds")

resultsNames(deseq_res)

#########################

# ALS

res_ALS_Zelic <- as.data.frame(results(deseq_res, name = "group_ALS_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_ALS_Zelic$t <- sign(res_ALS_Zelic$log2FoldChange) * -log10(res_ALS_Zelic$pvalue)
res_ALS_Zelic <- res_ALS_Zelic %>%
  arrange(desc(t))

ranking_ALS_Zelic <- res_ALS_Zelic$t
names(ranking_ALS_Zelic) <- res_ALS_Zelic$gene

gsea_ALS_Zelic <- fgsea(pathways = pathway, stats = ranking_ALS_Zelic, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_ALS_Zelic)

##################################################
##################################################
##################################################

# Kamath 2024

deseq_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Kamath_2022_dx_PD_deseq.rds")

resultsNames(deseq_res)

#########################

# PD

res_PD_Kamath <- as.data.frame(results(deseq_res, name = "group_PD.DLB_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_PD_Kamath$t <- sign(res_PD_Kamath$log2FoldChange) * -log10(res_PD_Kamath$pvalue)
res_PD_Kamath <- res_PD_Kamath %>%
  arrange(desc(t))

ranking_PD_Kamath <- res_PD_Kamath$t
names(ranking_PD_Kamath) <- res_PD_Kamath$gene

gsea_PD_Kamath <- fgsea(pathways = pathway, stats = ranking_PD_Kamath, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_PD_Kamath)

##################################################
##################################################
##################################################

# make overall summary

gsea_merged <- rbind(Green_res, Mathys_res, SEAAD_res, Rexach_res, Macnair_res, Petrescu_res, gsea_PD, 
                     Prudencio_res, gsea_PD_Martirosyan, gsea_ALS_Zelic, gsea_PD_Kamath)

gsea_merged$test <- c("Green_2024_PFC_CERAD", "Green_2024_PFC_Braak", "Green_2024_PFC_MCI", "Green_2024_PFC_AD_dementia",
                      "Mathys_2023_PFC_CERAD", "Mathys_2023_PFC_Braak", "Mathys_2023_PFC_MCI", "Mathys_2023_PFC_AD_dementia",
                      "Gabitto_2024_PFC_CERAD", "Gabitto_2024_PFC_Braak", "Gabitto_2024_PFC_AD_dementia",
                      "Rexach_2024_INS_AD", "Rexach_2024_INS_PSP", "Rexach_2024_INS_PiD",
                      "Macnair_2025_CTX_WM_NAWM", "Macnair_2025_CTX_WM_AL", "Macnair_2025_CTX_WM_CAL",
                      "Macnair_2025_CTX_WM_CIL", "Macnair_2025_CTX_WM_RL",
                      "Macnair_2025_CTX_GM_NAGM", "Macnair_2025_CTX_GM_GML",
                      "Petrescu_2025_BA44_ALS", "Petrescu_2025_BA44_ALS_borderline", "Petrescu_2025_BA44_ALSci",
                      "Petrescu_2025_BA46_ALS", "Petrescu_2025_BA46_ALS_borderline", "Petrescu_2025_BA46_ALSci",
                      "Wang_2024_SN_PD/PDD",
                      "Prudencio_2015_FC_sALS", "Prudencio_2015_FC_C9ALS",
                      "Prudencio_2015_CB_sALS", "Prudencio_2015_CB_C9ALS",
                      "Martirosyan_2024_SN_PD",
                      "Zelic_2025_SpC_ALS",
                      "Kamath_2022_SN_PD/DLB")

gsea_merged$padj <- p.adjust(gsea_merged$pval)

gsea_merged <- gsea_merged %>%
  dplyr::select(test, NES, size, pval, padj)

saveRDS(gsea_merged, file = "./analysis/microglia/DAM_signature_replication/DESeq_DAM_GSEA_summary.rds")

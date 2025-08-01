library(edgeR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(SeuratDisk)
library(EnvStats)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# prep counts data for PCA analysis

##################################################

# processed microglia object from Green et al., 2024, Nature (synapse ID = syn53693904)
Green_micro <- LoadH5Seurat("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/processed_data_from_study/microglia.h5Seurat")


# add pct.mito
Green_micro$pct.mito <- PercentageFeatureSet(Green_micro, pattern = "^MT")


# pseudobulk by donor
pb_cts <- as.data.frame(AggregateExpression(Green_micro, assays = "RNA", group.by = "individualID"))
pb_cts <- pb_cts[,-1]
colnames(pb_cts) <- gsub("RNA.", "", colnames(pb_cts))


# retain donors w/ at least 50 cells
t <- Green_micro@meta.data %>%
  group_by(individualID) %>%
  summarise(n = n())

donors_50cells <- t$individualID[t$n >= 50]

pb_cts <- pb_cts[,colnames(pb_cts) %in% donors_50cells]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(Green_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cpm.rds")
saveRDS(Green_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

ROSMAP_donor_meta <- read.csv("./metadata/ROSMAP_clinical.csv")
celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_celllevel_meta.rds")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cpm.rds")

##################################################

donor_meta <- ROSMAP_donor_meta[ROSMAP_donor_meta$individualID %in% colnames(pb_cpm), ]


celllevel_meta <- celllevel_meta[celllevel_meta$individualID %in% colnames(pb_cpm), ]
celllevel_meta <- celllevel_meta %>%
  distinct(individualID, .keep_all = T) %>%
  dplyr::select(individualID, batch)


donor_meta <- donor_meta %>%
  left_join(celllevel_meta, by = "individualID")


donor_meta <- donor_meta[c("individualID", "cogdx", "ceradsc", "braaksc", "age_death", "apoe_genotype", "pmi", "msex", "batch")]


donor_meta <- na.omit(donor_meta)


donor_meta <- donor_meta[donor_meta$cogdx %in% c("1", "2", "4"),]
donor_meta$cogdx <- ifelse(donor_meta$cogdx == "1", "HC",
                        ifelse(donor_meta$cogdx == "2", "MCI",
                               ifelse(donor_meta$cogdx == "4", "AD_dementia", NA)))
donor_meta$cogdx <- factor(donor_meta$cogdx, levels = c("HC", "MCI", "AD_dementia"))


donor_meta$cerad <- ifelse(donor_meta$ceradsc == "1", "Definite_AD",
                           ifelse(donor_meta$ceradsc == "2", "Probable_AD",
                                  ifelse(donor_meta$ceradsc == "3", "Possible_AD", "Not_AD")))
donor_meta$cerad <- factor(donor_meta$cerad, levels = c("Not_AD", "Possible_AD", "Probable_AD", "Definite_AD"))
donor_meta$cerad_numerical <- ifelse(donor_meta$ceradsc == "1", 1,
                                     ifelse(donor_meta$ceradsc == "2", 2/3,
                                            ifelse(donor_meta$ceradsc == "3", 1/3, 0)))


donor_meta$braak_scaled <- ifelse(donor_meta$braaksc == "0", 0,
                                  ifelse(donor_meta$braaksc == "1", 1/6,
                                         ifelse(donor_meta$braaksc == "2", 2/6,
                                                ifelse(donor_meta$braaksc == "3", 3/6,
                                                       ifelse(donor_meta$braaksc == "4", 4/6,
                                                              ifelse(donor_meta$braaksc == "5", 5/6, 1))))))
donor_meta$braak_scaled <- as.numeric(donor_meta$braak_scaled)


donor_meta$age_binned <- ifelse(donor_meta$age_death < 80, "<80", 
                                ifelse(donor_meta$age_death >= 80 & donor_meta$age_death < 90, "80-89", "90+"))
donor_meta$age_binned <- factor(donor_meta$age_binned, levels = c("<80", "80-89", "90+"))


donor_meta$apoe_genotype <- factor(donor_meta$apoe_genotype, levels = c("33", "22", "23", "24", "34", "44"))


donor_meta <- donor_meta %>%
  dplyr::select(individualID, age_death, age_binned, msex, pmi, batch, cogdx, apoe_genotype, braaksc, braak_scaled, cerad, cerad_numerical) %>%
  dplyr::rename(age = age_death, braak_score = braaksc)


# save
saveRDS(donor_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/donor_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
donor_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/donor_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])


# filter for usable donors
pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donor_meta$individualID,]


# run PCA
pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
pca_DAM_df <- as.data.frame(pca_DAM$x)


# join PCA df w/ metadata
pca_DAM_df <- pca_DAM_df %>%
  rownames_to_column(var = "individualID") %>%
  dplyr::select(individualID, PC1, PC2, PC3) %>%
  left_join(donor_meta, by = "individualID")


# test for outliers
rosnerTest(pca_DAM_df$PC1, k = 10)


# 7 outliers, remove & re-run PCA
donors_keep <- pca_DAM_df$individualID[pca_DAM_df$PC1 > -13]

pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donors_keep,]

pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
pca_DAM_df <- as.data.frame(pca_DAM$x)

pca_DAM_df <- pca_DAM_df %>%
  rownames_to_column(var = "individualID") %>%
  dplyr::select(individualID, PC1, PC1, PC2, PC3) %>%
  left_join(donor_meta, by = "individualID")

rosnerTest(pca_DAM_df$PC1, k = 10)
# no more outliers, proceed

#########################

sdev <- pca_DAM$sdev
variance <- sdev^2
variance_explained <- variance / sum(variance) * 100
# PC1: 17.85%

#########################

# save
saveRDS(pca_DAM_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final.rds")
saveRDS(variance_explained, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_pct_variance_explained.rds")

##################################################
##################################################
##################################################

# # run PCA analysis, just APOE3/3 and APOE3/4
# 
# pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cpm.rds")
# DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
# donor_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/donor_meta.rds")
# 
# ##################################################
# 
# # filter for DAM signature genes
# pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])
# 
# 
# donor_meta <- donor_meta[donor_meta$apoe_genotype %in% c("23", "33", "34") & donor_meta$cogdx == "AD_dementia", ]
# 
# 
# # filter for usable donors
# pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donor_meta$individualID,]
# 
# 
# # run PCA
# pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
# pca_DAM_df <- as.data.frame(pca_DAM$x)
# 
# 
# # join PCA df w/ metadata
# pca_DAM_df <- pca_DAM_df %>%
#   rownames_to_column(var = "individualID") %>%
#   dplyr::select(individualID, PC1, PC2, PC3) %>%
#   left_join(donor_meta, by = "individualID")
# 
# 
# # test for outliers
# rosnerTest(pca_DAM_df$PC1, k = 10)
# 
# 
# # 5 outliers, remove & re-run PCA
# donors_keep <- pca_DAM_df$individualID[pca_DAM_df$PC1 < 16.748703]
# 
# pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donors_keep,]
# 
# pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
# pca_DAM_df <- as.data.frame(pca_DAM$x)
# 
# pca_DAM_df <- pca_DAM_df %>%
#   rownames_to_column(var = "individualID") %>%
#   dplyr::select(individualID, PC1, PC1, PC2, PC3) %>%
#   left_join(donor_meta, by = "individualID")
# 
# rosnerTest(pca_DAM_df$PC1, k = 10)
# # no more outliers, proceed
# 
# #########################
# 
# sdev <- pca_DAM$sdev
# variance <- sdev^2
# variance_explained <- variance / sum(variance) * 100
# # PC1: 17.85%
# 
# #########################
# 
# # save
# saveRDS(pca_DAM_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final_APOE.rds")
# saveRDS(variance_explained, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_pct_variance_explained_APOE.rds")

##################################################
##################################################
##################################################

# run GLMs

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final.rds")
pca_DAM_df$pmi_scaled <- scale(pca_DAM_df$pmi)

##################################################

# CERAD

glm_cerad <- glm(formula = PC1 ~ cerad_numerical + msex + age_binned + pmi_scaled + batch, data = pca_DAM_df)
model_summary_cerad <- summary(glm_cerad)
coef_cerad <- model_summary_cerad$coefficients

#########################

# Braak

glm_braak <- glm(formula = PC1 ~ braak_scaled + msex + age_binned + pmi_scaled + batch, data = pca_DAM_df)
model_summary_braak <- summary(glm_braak)
coef_braak <- model_summary_braak$coefficients

#########################

# cognitive diagnosis

glm_cogdx <- glm(formula = PC1 ~ cogdx + msex + age_binned + pmi_scaled + batch, data = pca_DAM_df)
model_summary_cogdx <- summary(glm_cogdx)
coef_cogdx <- model_summary_cogdx$coefficients

#########################

# # APOE status
# 
# pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final_APOE.rds")
# pca_DAM_df$pmi_scaled <- scale(pca_DAM_df$pmi)
# 
# glm_apoe <- glm(formula = PC1 ~ apoe_genotype + msex + age_binned + pmi_scaled + batch, data = pca_DAM_df)
# model_summary_apoe <- summary(glm_apoe)
# coef_apoe <- model_summary_apoe$coefficients

#########################

# merge

glm_summary_merged <- as.data.frame(rbind(coef_cerad, coef_braak, coef_cogdx))

glm_summary_merged <- glm_summary_merged %>%
  rownames_to_column(var = "test")
glm_summary_merged <- glm_summary_merged[grepl("braak|cerad|dx", glm_summary_merged$test),]
glm_summary_merged$dataset <- "Green_2024"

rownames(glm_summary_merged) <- NULL

#########################

# save
saveRDS(glm_summary_merged, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_GLM_summary.rds")

##################################################
##################################################
##################################################

# generate data to use w/ variancePartition

varpart_gex <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_pb_cpm.rds")
varpart_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/donor_meta.rds")
pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final.rds")

#########################

rownames(varpart_meta) <- NULL
varpart_meta <- varpart_meta %>%
  column_to_rownames(var = "individualID")
varpart_meta <- varpart_meta[rownames(varpart_meta) %in% pca_DAM_df$individualID, ]


varpart_gex <- varpart_gex[, colnames(varpart_gex) %in% rownames(varpart_meta)]
varpart_gex <- varpart_gex[, rownames(varpart_meta)]


# save
saveRDS(varpart_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/varpart_meta.rds")
saveRDS(varpart_gex, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/varpart_GEx.rds")

##################################################
##################################################
##################################################

# prep filtered single cell counts matrix, cell-level metadata, and donor metadata for further use

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final.rds")
celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_celllevel_meta.rds")
sc_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_sc_cts.rds")

#########################

sc_cts <- sc_cts[, colnames(sc_cts) %in% rownames(celllevel_meta)[celllevel_meta$individualID %in% pca_DAM_df$individualID]]


celllevel_meta <- celllevel_meta[celllevel_meta$individualID %in% pca_DAM_df$individualID, ]
celllevel_meta <- celllevel_meta %>%
  dplyr::select(individualID, nCount_RNA, pct.mito) %>%
  rownames_to_column(var = "barcode") %>%
  left_join(pca_DAM_df, by = "individualID") %>%
  column_to_rownames(var = "barcode")
  

sc_cts <- sc_cts[, rownames(celllevel_meta)]


identical(rownames(celllevel_meta), colnames(sc_cts))

#########################

# save
saveRDS(sc_cts, "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_sc_cts_FILTERED.rds")
saveRDS(celllevel_meta, "./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/microglia_celllevel_meta_FILTERED.rds")

library(edgeR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(SeuratDisk)
library(EnvStats)
library(schard)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# prep counts data for PCA analysis

##################################################

# microglia from SEA-AD processed PFC object, downloaded from AWS 
# https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad
# subset out "Micro-PVM" subclass with scanpy
SEAAD_micro <- schard::h5ad2seurat("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/processed_data_from_study/SEAAD_PFC_microglia.h5ad")


# add nCount_RNA
SEAAD_micro$nCount_RNA <- colSums(SEAAD_micro@assays$RNA$counts)
SEAAD_micro$pct.mito <- PercentageFeatureSet(SEAAD_micro, pattern = "^MT")


# pseudobulk by donor
pb_cts <- as.data.frame(AggregateExpression(SEAAD_micro, assays = "RNA", group.by = "Donor ID"))
colnames(pb_cts) <- gsub("RNA.", "", colnames(pb_cts))


# retain donors w/ at least 50 cells
t <- SEAAD_micro@meta.data %>%
  group_by(`Donor ID`) %>%
  summarise(n = n())

donors_50cells <- t$`Donor ID`[t$n >= 50]

pb_cts <- pb_cts[,colnames(pb_cts) %in% donors_50cells]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(SEAAD_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_pb_cpm.rds")
saveRDS(SEAAD_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

SEAAD_donor_meta <- readRDS("/data/ADRD/SEA_AD/metadata/donor_meta_final.rds")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_pb_cpm.rds")

##################################################

donor_meta <- SEAAD_donor_meta[SEAAD_donor_meta$donor %in% colnames(pb_cpm), ]


donor_meta <- donor_meta[donor_meta$severely_affected_donor == "N", ]


donor_meta <- donor_meta[c("donor", "age_at_death", "sex", "APOE_genotype", "braak_stage", "cerad_numerical", "AD_dementia_status", "pmi")]


donor_meta$cogdx <- ifelse(donor_meta$AD_dementia_status == "dementia", "AD_dementia", "HC")
donor_meta$cogdx <- factor(donor_meta$cogdx, levels = c("HC", "AD_dementia"))


donor_meta$cerad <- ifelse(donor_meta$cerad_numerical == "3", "Definite_AD",
                           ifelse(donor_meta$cerad_numerical == "2", "Probable_AD",
                                  ifelse(donor_meta$cerad_numerical == "1", "Possible_AD", "No_AD")))
donor_meta$cerad <- factor(donor_meta$cerad, levels = c("No_AD", "Possible_AD", "Probable_AD", "Definite_AD"))
donor_meta$cerad_numerical <- ifelse(donor_meta$cerad_numerical == "3", 1,
                                     ifelse(donor_meta$cerad_numerical == "2", 2/3,
                                            ifelse(donor_meta$cerad_numerical == "1", 1/3, 0)))



donor_meta$braak_scaled <- ifelse(donor_meta$braak_stage == "0", 0,
                                  ifelse(donor_meta$braak_stage == "1", 1/6,
                                         ifelse(donor_meta$braak_stage == "2", 2/6,
                                                ifelse(donor_meta$braak_stage == "3", 3/6,
                                                       ifelse(donor_meta$braak_stage == "4", 4/6,
                                                              ifelse(donor_meta$braak_stage == "5", 5/6, 1))))))
donor_meta$braak_scaled <- as.numeric(donor_meta$braak_scaled)


donor_meta$APOE_genotype <- factor(donor_meta$APOE_genotype, levels = c("33", "22", "23", "24", "34", "44"))


donor_meta <- donor_meta %>%
  dplyr::select(donor, age_at_death, sex, pmi, cogdx, APOE_genotype, braak_stage, braak_scaled, cerad, cerad_numerical) %>%
  dplyr::rename(age = age_at_death, braak_score = braak_stage, apoe_genotype = APOE_genotype)


# save
saveRDS(donor_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/donor_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
donor_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/donor_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])


# filter for usable donors
pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donor_meta$donor,]


# run PCA
pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
pca_DAM_df <- as.data.frame(pca_DAM$x)


# join PCA df w/ metadata
pca_DAM_df <- pca_DAM_df %>%
  rownames_to_column(var = "donor") %>%
  dplyr::select(donor, PC1, PC2, PC3) %>%
  left_join(donor_meta, by = "donor")


# test for outliers
rosnerTest(pca_DAM_df$PC1, k = 10)


# 2 outliers, remove & re-run PCA
donors_keep <- pca_DAM_df$donor[pca_DAM_df$PC1 > -16]

pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donors_keep,]


pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
pca_DAM_df <- as.data.frame(pca_DAM$x)

pca_DAM_df <- pca_DAM_df %>%
  rownames_to_column(var = "donor") %>%
  dplyr::select(donor, PC1, PC1, PC2, PC3) %>%
  left_join(donor_meta, by = "donor")

rosnerTest(pca_DAM_df$PC1, k = 10)
# no more outliers, proceed

#########################

sdev <- pca_DAM$sdev
variance <- sdev^2
variance_explained <- variance / sum(variance) * 100
# PC1: 23.13%

#########################

# save
saveRDS(pca_DAM_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_DAM_final.rds")
saveRDS(variance_explained, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_pct_variance_explained.rds")

##################################################
##################################################
##################################################

# run GLMs

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_DAM_final.rds")
pca_DAM_df$pmi_scaled <- scale(pca_DAM_df$pmi)
pca_DAM_df$age_scaled <- scale(pca_DAM_df$age)

##################################################

# CERAD

glm_cerad <- glm(formula = PC1 ~ cerad_numerical + sex + age_scaled + pmi_scaled, data = pca_DAM_df)
model_summary_cerad <- summary(glm_cerad)
coef_cerad <- model_summary_cerad$coefficients

#########################

# Braak

glm_braak <- glm(formula = PC1 ~ braak_scaled + sex + age_scaled + pmi_scaled, data = pca_DAM_df)
model_summary_braak <- summary(glm_braak)
coef_braak <- model_summary_braak$coefficients

#########################

# cognitive diagnosis

glm_cogdx <- glm(formula = PC1 ~ cogdx + sex + age_scaled + pmi_scaled, data = pca_DAM_df)
model_summary_cogdx <- summary(glm_cogdx)
coef_cogdx <- model_summary_cogdx$coefficients

#########################

# merge

glm_summary_merged <- as.data.frame(rbind(coef_cerad, coef_braak, coef_cogdx))

glm_summary_merged <- glm_summary_merged %>%
  rownames_to_column(var = "test")
glm_summary_merged <- glm_summary_merged[grepl("braak|cerad|dx", glm_summary_merged$test),]
glm_summary_merged$dataset <- "SEAAD_PFC"

rownames(glm_summary_merged) <- NULL

#########################

# save
saveRDS(glm_summary_merged, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_GLM_summary.rds")

##################################################
##################################################
##################################################

# generate data to use w/ variancePartition

varpart_gex <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_pb_cpm.rds")
varpart_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/donor_meta.rds")
pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_DAM_final.rds")

#########################

rownames(varpart_meta) <- NULL
varpart_meta <- varpart_meta %>%
  column_to_rownames(var = "donor")
varpart_meta <- varpart_meta[rownames(varpart_meta) %in% pca_DAM_df$donor, ]


varpart_gex <- varpart_gex[, colnames(varpart_gex) %in% rownames(varpart_meta)]
varpart_gex <- varpart_gex[, rownames(varpart_meta)]


# save
saveRDS(varpart_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/varpart_meta.rds")
saveRDS(varpart_gex, file = "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/varpart_GEx.rds")

##################################################
##################################################
##################################################

# prep filtered single cell counts matrix, cell-level metadata, and donor metadata for further use

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/PCA_DAM_final.rds")
celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_celllevel_meta.rds")
sc_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_sc_cts.rds")

#########################

sc_cts <- sc_cts[, colnames(sc_cts) %in% rownames(celllevel_meta)[celllevel_meta$`Donor ID` %in% pca_DAM_df$donor]]


celllevel_meta <- celllevel_meta[celllevel_meta$`Donor ID` %in% pca_DAM_df$donor, ]
celllevel_meta$`Donor ID` <- as.character(celllevel_meta$`Donor ID`)
celllevel_meta <- celllevel_meta %>%
  dplyr::rename(donor = `Donor ID`) %>%
  dplyr::select(donor, nCount_RNA, pct.mito) %>%
  rownames_to_column(var = "barcode") %>%
  left_join(pca_DAM_df, by = "donor") %>%
  column_to_rownames(var = "barcode")


sc_cts <- sc_cts[, rownames(celllevel_meta)]


identical(rownames(celllevel_meta), colnames(sc_cts))

#########################

# save
saveRDS(sc_cts, "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_sc_cts_FILTERED.rds")
saveRDS(celllevel_meta, "./analysis/microglia/DAM_signature_replication/individual_studies/SEAAD_PFC/microglia_celllevel_meta_FILTERED.rds")

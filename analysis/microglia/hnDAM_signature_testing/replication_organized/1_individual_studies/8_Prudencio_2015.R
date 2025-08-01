library(edgeR)
library(tidyverse)
library(ggplot2)
library(EnvStats)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# prep counts data for PCA analysis

##################################################

# note this is bulk RNAseq!!!

# counts matrix downloaded from GSE67196, donor meta from Supp Table 1 of paper
cts <- read.delim("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/processed_data_from_study/GSE67196_Petrucelli2015_ALS_genes.rawcount.txt")


cts <- cts[, -1]
cts$GeneID <- make.unique(cts$GeneID)
cts <- cts %>%
  column_to_rownames(var = "GeneID")


# split counts fcx/cereb
cts_fcx <- cts[, grepl("fcx", colnames(cts))]
cts_cereb <- cts[, grepl("cereb", colnames(cts))]


colnames(cts_cereb) <- gsub("ALS", "", colnames(cts_cereb))
colnames(cts_cereb) <- sub("^0*", "", colnames(cts_cereb))
colnames(cts_cereb) <- gsub("_cereb", "", colnames(cts_cereb))

colnames(cts_fcx) <- gsub("ALS", "", colnames(cts_fcx))
colnames(cts_fcx) <- sub("^0*", "", colnames(cts_fcx))
colnames(cts_fcx) <- gsub("_fcx", "", colnames(cts_fcx))


# TMM-normalize counts
dge_fcx <- DGEList(cts_fcx)
dge_fcx <- calcNormFactors(dge_fcx, method = "TMM")
cpm_fcx <- as.data.frame(edgeR::cpm(dge_fcx))

dge_cereb <- DGEList(cts_cereb)
dge_cereb <- calcNormFactors(dge_cereb, method = "TMM")
cpm_cereb <- as.data.frame(edgeR::cpm(dge_cereb))


# save
saveRDS(cts_fcx, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/cts_fcx.rds")
saveRDS(cts_cereb, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/cts_cereb.rds")
saveRDS(cpm_fcx, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/cpm_fcx.rds")
saveRDS(cpm_cereb, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/cpm_cereb.rds")

##################################################
##################################################
##################################################

# run PCA analysis

cpm_fcx <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/cpm_fcx.rds")
cpm_cereb <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/cpm_cereb.rds")

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
donor_meta <- read.csv("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/Prudencio_2015_donor_meta.csv")
donor_meta$donor <- as.character(donor_meta$donor)

##################################################

# filter for DAM signature genes
cpm_fcx_DAM <- t(cpm_fcx[rownames(cpm_fcx) %in% c(DAM_genes), ])
cpm_cereb_DAM <- t(cpm_cereb[rownames(cpm_cereb) %in% c(DAM_genes), ])

##################################################

# run PCA

#########################

# frontal cortex

pca_DAM_fcx <- prcomp(cpm_fcx_DAM, center = TRUE, scale. = TRUE)
pca_DAM_fcx_df <- as.data.frame(pca_DAM_fcx$x)


# join PCA df w/ metadata
pca_DAM_fcx_df <- pca_DAM_fcx_df %>%
  rownames_to_column(var = "donor") %>%
  dplyr::select(donor, PC1, PC2, PC3) %>%
  left_join(donor_meta, by = "donor")


# test for outliers
rosnerTest(pca_DAM_fcx_df$PC1, k = 10)
# no outliers, proceed


pca_DAM_fcx_df$age_scaled <- scale(pca_DAM_fcx_df$age_at_death)
pca_DAM_fcx_df$pmi_scaled <- scale(pca_DAM_fcx_df$pmi)
pca_DAM_fcx_df$group <- factor(pca_DAM_fcx_df$group, levels = c("HC", "sALS", "C9ALS"))

#########################

sdev_fcx <- pca_DAM_fcx$sdev
variance_fcx <- sdev_fcx^2
variance_explained_fcx <- variance_fcx / sum(variance_fcx) * 100
# PC1: 32.93%

#########################

# save
saveRDS(pca_DAM_fcx_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_DAM_final_fcx.rds")
saveRDS(variance_explained_fcx, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_pct_variance_explained_fcx.rds")

##################################################

# cerebellum

pca_DAM_cereb <- prcomp(cpm_cereb_DAM, center = TRUE, scale. = TRUE)
pca_DAM_cereb_df <- as.data.frame(pca_DAM_cereb$x)


# join PCA df w/ metadata
pca_DAM_cereb_df <- pca_DAM_cereb_df %>%
  rownames_to_column(var = "donor") %>%
  dplyr::select(donor, PC1, PC2, PC3) %>%
  left_join(donor_meta, by = "donor")


# test for outliers
rosnerTest(pca_DAM_cereb_df$PC1, k = 10)
# no outliers, proceed


pca_DAM_cereb_df$age_scaled <- scale(pca_DAM_cereb_df$age_at_death)
pca_DAM_cereb_df$pmi_scaled <- scale(pca_DAM_cereb_df$pmi)
pca_DAM_cereb_df$group <- factor(pca_DAM_cereb_df$group, levels = c("HC", "sALS", "C9ALS"))

#########################

sdev_cereb <- pca_DAM_cereb$sdev
variance_cereb <- sdev_cereb^2
variance_explained_cereb <- variance_cereb / sum(variance_cereb) * 100
# PC1: 26.00%

#########################

# save
saveRDS(pca_DAM_cereb_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_DAM_final_cereb.rds")
saveRDS(variance_explained_cereb, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_pct_variance_explained_cereb.rds")

##################################################
##################################################
##################################################

# run GLMs

pca_DAM_fcx_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_DAM_final_fcx.rds")
pca_DAM_cereb_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_DAM_final_cereb.rds")

##################################################

glm_fcx <- glm(formula = PC1 ~ group + sex + age_scaled + pmi_scaled, data = pca_DAM_fcx_df)
model_summary_fcx <- summary(glm_fcx)
coef_fcx <- model_summary_fcx$coefficients
coef_fcx <- coef_fcx %>%
  as.data.frame(.) %>%
  rownames_to_column(var = "test")
coef_fcx$test <- paste0(coef_fcx$test, "_fcx")


glm_cereb <- glm(formula = PC1 ~ group + sex + age_scaled + pmi_scaled, data = pca_DAM_cereb_df)
model_summary_cereb <- summary(glm_cereb)
coef_cereb <- model_summary_cereb$coefficients
coef_cereb <- coef_cereb %>%
  as.data.frame(.) %>%
  rownames_to_column(var = "test")
coef_cereb$test <- paste0(coef_cereb$test, "_cereb")

#########################

# merge

glm_summary_merged <- as.data.frame(rbind(coef_fcx, coef_cereb))

glm_summary_merged <- glm_summary_merged[grepl("group", glm_summary_merged$test),]
glm_summary_merged$dataset <- "Prudencio_2015"

rownames(glm_summary_merged) <- NULL

#########################

# save
saveRDS(glm_summary_merged, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Prudencio_2015/PCA_GLM_summary.rds")

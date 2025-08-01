library(edgeR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(SeuratDisk)
library(EnvStats)
library(Matrix)
library(lmerTest)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# prep counts data for PCA analysis

##################################################

# processed data downloaded from zenodo: https://zenodo.org/records/8338963
cts <- readMM("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/processed_data_from_study/ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx")
cts <- as(cts, "CsparseMatrix")

genes <- read.delim("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/processed_data_from_study/ms_lesions_snRNAseq_row_data_2023-09-12.txt", sep = ",")
barcodes <- read.delim("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/processed_data_from_study/ms_lesions_snRNAseq_col_data_2023-09-12.txt", sep = ",")

colnames(cts) <- barcodes$cell_id
rownames(cts) <- make.unique(genes$symbol)

barcodes <- barcodes[barcodes$type_broad == "Microglia", ]
barcodes <- barcodes[barcodes$exclude_pseudobulk == "FALSE", ]

cts <- cts[, colnames(cts) %in% barcodes$cell_id]

rownames(barcodes) <- NULL
celllevel_meta <- barcodes %>%
  column_to_rownames(var = "cell_id")

Macnair_micro <- CreateSeuratObject(counts = cts, meta.data = celllevel_meta)

Macnair_micro$pct.mito <- PercentageFeatureSet(Macnair_micro, pattern = "^MT", assay = "RNA")

#########################

# pseudobulk by sample
pb_cts <- as.data.frame(AggregateExpression(Macnair_micro, assays = "RNA", group.by = "sample_id_anon"))
colnames(pb_cts) <- gsub("RNA.", "", colnames(pb_cts))
colnames(pb_cts) <- gsub("\\.", "-", colnames(pb_cts))


# retain samples w/ at least 50 cells
t <- Macnair_micro@meta.data %>%
  group_by(sample_id_anon) %>%
  summarise(n = n())

donors_50cells <- t$sample_id_anon[t$n >= 50]

pb_cts <- pb_cts[,colnames(pb_cts) %in% donors_50cells]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(Macnair_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_pb_cpm.rds")
saveRDS(Macnair_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_celllevel_meta.rds")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_pb_cpm.rds")

##################################################

sample_meta <- celllevel_meta %>%
  distinct(sample_id_anon, .keep_all = T)


sample_meta <- sample_meta[sample_meta$sample_id_anon %in% colnames(pb_cpm), ]


sample_meta <- sample_meta %>%
  dplyr::select(sample_id_anon, individual_id_anon, lesion_type, seq_pool, matter, sample_source, diagnosis, sex, age_at_death, pmi_minutes)


sample_meta$pmi <- sample_meta$pmi_minutes / 60


sample_meta <- sample_meta %>%
  dplyr::rename(age = age_at_death) %>%
  dplyr::select(-pmi_minutes)


# save
saveRDS(sample_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/sample_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
sample_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/sample_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])

##################################################

# white matter samples

pb_cpm_DAM_WM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% sample_meta$sample_id_anon[sample_meta$matter == "WM"], ]


# run PCA
pca_DAM_WM <- prcomp(pb_cpm_DAM_WM, center = TRUE, scale. = TRUE)
pca_DAM_df_WM <- as.data.frame(pca_DAM_WM$x)


rosnerTest(pca_DAM_df_WM$PC1, k = 10)
# no outliers, continue


pca_DAM_df_WM <- pca_DAM_df_WM %>%
  rownames_to_column(var = "sample_id_anon") %>%
  dplyr::select(sample_id_anon, PC1, PC2, PC3) %>%
  left_join(sample_meta, by = "sample_id_anon")


pca_DAM_df_WM$age_scaled <- scale(pca_DAM_df_WM$age)
pca_DAM_df_WM$pmi_scaled <- scale(pca_DAM_df_WM$pmi)
pca_DAM_df_WM$lesion_type <- factor(pca_DAM_df_WM$lesion_type, levels = c("WM", "NAWM", "AL", "CAL", "CIL", "RL"))

#########################

sdev_WM <- pca_DAM_WM$sdev
variance_WM <- sdev_WM^2
variance_explained_WM <- variance_WM / sum(variance_WM) * 100
# PC1: 30.64%

#########################

# save
saveRDS(pca_DAM_df_WM, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_DAM_final_WM.rds")
saveRDS(variance_explained_WM, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_pct_variance_explained_WM.rds")

##################################################

# grey matter samples

pb_cpm_DAM_GM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% sample_meta$sample_id_anon[sample_meta$matter == "GM"], ]
pb_cpm_DAM_GM <- pb_cpm_DAM_GM[, !colnames(pb_cpm_DAM_GM) == "CTSD"]

# run PCA
pca_DAM_GM <- prcomp(pb_cpm_DAM_GM, center = TRUE, scale. = TRUE)
pca_DAM_df_GM <- as.data.frame(pca_DAM_GM$x)


pca_DAM_df_GM <- pca_DAM_df_GM %>%
  rownames_to_column(var = "sample_id_anon") %>%
  dplyr::select(sample_id_anon, PC1, PC2, PC3) %>%
  left_join(sample_meta, by = "sample_id_anon")


rosnerTest(pca_DAM_df_GM$PC1, k = 10)


# 1 outlier, remove & re-run PCA
samples_keep <- pca_DAM_df_GM$sample_id_anon[pca_DAM_df_GM$PC1 > -20]

pb_cpm_DAM_GM <- pb_cpm_DAM_GM[rownames(pb_cpm_DAM_GM) %in% samples_keep,]

pca_DAM_GM <- prcomp(pb_cpm_DAM_GM, center = TRUE, scale. = TRUE)
pca_DAM_df_GM <- as.data.frame(pca_DAM_GM$x)

pca_DAM_df_GM <- pca_DAM_df_GM %>%
  rownames_to_column(var = "sample_id_anon") %>%
  dplyr::select(sample_id_anon, PC1, PC2, PC3) %>%
  left_join(sample_meta, by = "sample_id_anon")

rosnerTest(pca_DAM_df_GM$PC1, k = 10)
# no more outliers, continue


pca_DAM_df_GM$age_scaled <- scale(pca_DAM_df_GM$age)
pca_DAM_df_GM$pmi_scaled <- scale(pca_DAM_df_GM$pmi)
pca_DAM_df_GM$lesion_type <- factor(pca_DAM_df_GM$lesion_type, levels = c("GM", "NAGM", "GML"))

#########################

sdev_GM <- pca_DAM_GM$sdev
variance_GM <- sdev_GM^2
variance_explained_GM <- variance_GM / sum(variance_GM) * 100
# PC1: 25.3%

#########################

# save
saveRDS(pca_DAM_df_GM, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_DAM_final_GM.rds")
saveRDS(variance_explained_GM, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_pct_variance_explained_GM.rds")

##################################################
##################################################
##################################################

# run GLMs

pca_DAM_df_WM <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_DAM_final_WM.rds")
pca_DAM_df_GM <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_DAM_final_GM.rds")

##################################################

glm_WM <- lmer(formula = PC1 ~ lesion_type + sex + age_scaled + pmi_scaled + (1|individual_id_anon), data = pca_DAM_df_WM)
model_summary_WM <- summary(glm_WM)
coef_WM <- model_summary_WM$coefficients

glm_GM <- lmer(formula = PC1 ~ lesion_type + sex + age_scaled + pmi_scaled + (1|individual_id_anon), data = pca_DAM_df_GM)
model_summary_GM <- summary(glm_GM)
coef_GM <- model_summary_GM$coefficients

#########################

# merge

glm_summary_merged <- as.data.frame(rbind(coef_WM, coef_GM))

glm_summary_merged <- glm_summary_merged %>%
  rownames_to_column(var = "test")
glm_summary_merged <- glm_summary_merged[grepl("lesion", glm_summary_merged$test),]
glm_summary_merged$dataset <- "Macnair_2025"

rownames(glm_summary_merged) <- NULL

#########################

# save
saveRDS(glm_summary_merged, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_GLM_summary.rds")

##################################################
##################################################
##################################################

# prep filtered single cell counts matrix, cell-level metadata, and donor metadata for further use

sample_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/sample_meta.rds")
celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_celllevel_meta.rds")
sc_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_sc_cts.rds")

#########################

sc_cts <- sc_cts[, colnames(sc_cts) %in% rownames(celllevel_meta)[celllevel_meta$sample_id_anon %in% sample_meta$sample_id_anon]]


celllevel_meta <- celllevel_meta[celllevel_meta$sample_id_anon %in% sample_meta$sample_id_anon, ]
celllevel_meta <- celllevel_meta %>%
  dplyr::select(sample_id_anon, nCount_RNA, pct.mito) %>%
  rownames_to_column(var = "barcode") %>%
  left_join(sample_meta, by = "sample_id_anon") %>%
  column_to_rownames(var = "barcode")


sc_cts <- sc_cts[, rownames(celllevel_meta)]


identical(rownames(celllevel_meta), colnames(sc_cts))

#########################

# save
saveRDS(sc_cts, "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_sc_cts_FILTERED.rds")
saveRDS(celllevel_meta, "./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/microglia_celllevel_meta_FILTERED.rds")

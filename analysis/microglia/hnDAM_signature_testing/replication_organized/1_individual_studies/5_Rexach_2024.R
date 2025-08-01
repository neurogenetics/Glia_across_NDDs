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

# processed data and metadata downloaded from synapse (syn52082747)
seurat <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/processed_data_from_study/pci_seurat.rds")
meta <- read.csv("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/processed_data_from_study/liger_subcluster_metadata_v2.csv")


# filter for microglia & donors w/ 50+ microglia, and cells from insula
meta <- meta[meta$liger_clusters_stringent_retained == T, ]
meta <- meta[meta$cluster_cell_type == "microglia", ]
meta <- meta[meta$region == "insula", ]

t <- meta %>%
  group_by(library_id) %>%
  summarise(n = n())
donors_keep <- t$library_id[t$n >= 50]

meta <- meta[meta$library_id %in% donors_keep, ]

donors_keep <- gsub("_", ".", donors_keep)

cts <- seurat@assays$RNA$counts

micro_cts <- cts[, colnames(cts) %in% meta$UMI]

rownames(meta) <- NULL
meta <- meta %>%
  column_to_rownames(var = "UMI")

Rexach_micro <- CreateSeuratObject(counts = micro_cts, meta.data = meta)


# pseudobulk by donor
pb_cts <- as.data.frame(AggregateExpression(Rexach_micro, assays = "RNA", group.by = "library_id"))
colnames(pb_cts) <- gsub("RNA.", "", colnames(pb_cts))

pb_cts <- pb_cts[, colnames(pb_cts) %in% donors_keep]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(Rexach_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_pb_cpm.rds")
saveRDS(Rexach_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_celllevel_meta.rds")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_pb_cpm.rds")

##################################################

sample_meta <- celllevel_meta %>%
  distinct(library_id, .keep_all = T)


sample_meta <- sample_meta %>%
  dplyr::select(library_id, prep, autopsy_id, finalsite, npdx1, age, sex, pmi)


rownames(sample_meta) <- NULL


sample_meta$npdx1 <- ifelse(sample_meta$npdx1 %in% c("Control", "Normal"), "HC", sample_meta$npdx1)
sample_meta$npdx1 <- factor(sample_meta$npdx1, levels = c("HC", "AD", "PSP", "PiD"))


# save
saveRDS(sample_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/sample_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
sample_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/sample_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])

pb_cpm_DAM <- pb_cpm_DAM[, !colnames(pb_cpm_DAM) == "CTSD"]


# run PCA
pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
pca_DAM_df <- as.data.frame(pca_DAM$x)

rownames(pca_DAM_df) <- gsub("\\.", "_", rownames(pca_DAM_df))

# join PCA df w/ metadata
pca_DAM_df <- pca_DAM_df %>%
  rownames_to_column(var = "library_id") %>%
  dplyr::select(library_id, PC1, PC2, PC3) %>%
  left_join(sample_meta, by = "library_id")


rosnerTest(pca_DAM_df$PC1, k = 10)
# no outliers, continue


pca_DAM_df$npdx1 <- factor(pca_DAM_df$npdx1, levels = c("HC", "AD", "PSP", "PiD"))

#########################

sdev <- pca_DAM$sdev
variance <- sdev^2
variance_explained <- variance / sum(variance) * 100
# PC1: 22.43%

#########################

# save
saveRDS(pca_DAM_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/PCA_DAM_final.rds")
saveRDS(variance_explained, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/PCA_pct_variance_explained.rds")

##################################################
##################################################
##################################################

# run GLMs

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/PCA_DAM_final.rds")
pca_DAM_df$pmi_scaled <- scale(pca_DAM_df$pmi)
pca_DAM_df$age_scaled <- scale(pca_DAM_df$age)

##################################################

glm_npdx <- glm(formula = PC1 ~ npdx1 + sex + age_scaled + pmi_scaled + finalsite + prep, data = pca_DAM_df)
model_summary_npdx <- summary(glm_npdx)
coef_npdx <- model_summary_npdx$coefficients

#########################

glm_summary <- coef_npdx %>%
  as.data.frame(.) %>%
  rownames_to_column(var = "test")
glm_summary <- glm_summary[grepl("npdx1", glm_summary$test),]
glm_summary$dataset <- "Rexach_2024"

rownames(glm_summary) <- NULL

#########################

# save
saveRDS(glm_summary, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/PCA_GLM_summary.rds")

##################################################
##################################################
##################################################

# prep filtered single cell counts matrix, cell-level metadata, and donor metadata for further use

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/PCA_DAM_final.rds")
celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_celllevel_meta.rds")
sc_cts <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_sc_cts.rds")

#########################

sc_cts <- sc_cts[, colnames(sc_cts) %in% rownames(celllevel_meta)[celllevel_meta$library_id %in% pca_DAM_df$library_id]]


celllevel_meta <- celllevel_meta[celllevel_meta$library_id %in% pca_DAM_df$library_id, ]
celllevel_meta <- celllevel_meta %>%
  dplyr::select(library_id, nCount_RNA, percent.mt) %>%
  dplyr::rename(pct.mito = percent.mt) %>%
  rownames_to_column(var = "barcode") %>%
  left_join(pca_DAM_df, by = "library_id") %>%
  column_to_rownames(var = "barcode")


sc_cts <- sc_cts[, rownames(celllevel_meta)]


identical(rownames(celllevel_meta), colnames(sc_cts))

#########################

# save
saveRDS(sc_cts, "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_sc_cts_FILTERED.rds")
saveRDS(celllevel_meta, "./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/microglia_celllevel_meta_FILTERED.rds")

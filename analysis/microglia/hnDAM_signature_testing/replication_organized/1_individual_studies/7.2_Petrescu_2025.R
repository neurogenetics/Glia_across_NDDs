library(edgeR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(EnvStats)
library(schard)
library(AnnotationHub)
library(ensembldb)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# prep counts data for PCA analysis

##################################################

# processed microglia data from Petrescu et al., 2025 (preprint) -- data from GEO (GSE290359_BA44_BA46_COUNTS.h5ad)
# filtered microglia from main object w/ scanpy
seurat <- schard::h5ad2seurat("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_subset.h5ad")


# get counts -- need to replace ensembl IDs w/ gene names
cts <- LayerData(seurat, layer = "counts", assay = "RNA")

cellranger_map <- read.csv("/data/ADRD/human_brain_atlasing/3_fastq_processing/cellranger_features_2020.tsv", sep = "\t", header = F)
cellranger_map$V2 <- make.unique(cellranger_map$V2)
gene_mapping <- setNames(cellranger_map$V2, cellranger_map$V1)

rownames(cts) <- unname(gene_mapping[rownames(cts)])


Petrescu_micro <- CreateSeuratObject(counts = cts, meta.data = seurat@meta.data)


# add pct.mito
Petrescu_micro$pct.mito <- PercentageFeatureSet(Petrescu_micro, pattern = "^MT")


# pseudobulk by sample
pb_cts <- as.data.frame(AggregateExpression(Petrescu_micro, assays = "RNA", group.by = "library_name"))
colnames(pb_cts) <- gsub("RNA.", "", colnames(pb_cts))


# retain donors w/ at least 50 cells
Petrescu_micro$library_name <- gsub("-", ".", Petrescu_micro$library_name)
Petrescu_micro$library_name <- gsub("_", ".", Petrescu_micro$library_name)

t <- Petrescu_micro@meta.data %>%
  group_by(library_name) %>%
  summarise(n = n())

donors_50cells <- t$library_name[t$n >= 50]

pb_cts <- pb_cts[,colnames(pb_cts) %in% donors_50cells]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(Petrescu_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_pb_cpm.rds")
saveRDS(Petrescu_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_celllevel_meta.rds")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_pb_cpm.rds")

##################################################

sample_meta <- celllevel_meta %>%
  distinct(library_name, .keep_all = T) %>%
  dplyr::select(library_name, Classifcation, Sex, Age, Tissue)

rownames(sample_meta) <- NULL

sample_meta <- sample_meta[sample_meta$library_name %in% colnames(pb_cpm), ]

sample_meta$Classifcation <- factor(sample_meta$Classifcation, levels = c("Control", "ALS", "ALS - borderline", "ALSci"))

# save
saveRDS(sample_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/sample_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
sample_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/sample_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])

##################################################

# BA44

pb_cpm_DAM_BA44 <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% sample_meta$library_name[sample_meta$Tissue == "BA44"], ]

# run PCA
pca_DAM_BA44 <- prcomp(pb_cpm_DAM_BA44, center = TRUE, scale. = TRUE)
pca_DAM_df_BA44 <- as.data.frame(pca_DAM_BA44$x)


rosnerTest(pca_DAM_df_BA44$PC1, k = 5)
# no outliers, continue

pca_DAM_df_BA44 <- pca_DAM_df_BA44 %>%
  rownames_to_column(var = "library_name") %>%
  dplyr::select(library_name, PC1, PC2, PC3) %>%
  left_join(sample_meta, by = "library_name")


pca_DAM_df_BA44$Classifcation <- factor(pca_DAM_df_BA44$Classifcation, levels = c("Control", "ALS", "ALS - borderline", "ALSci"))

#########################

sdev_BA44 <- pca_DAM_BA44$sdev
variance_BA44 <- sdev_BA44^2
variance_explained_BA44 <- variance_BA44 / sum(variance_BA44) * 100
# PC1: 26.08%

#########################

# save
saveRDS(pca_DAM_df_BA44, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_DAM_final_BA44.rds")
saveRDS(variance_explained_BA44, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_pct_variance_explained_BA44.rds")

##################################################

# BA46

pb_cpm_DAM_BA46 <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% sample_meta$library_name[sample_meta$Tissue == "BA46"], ]

# run PCA
pca_DAM_BA46 <- prcomp(pb_cpm_DAM_BA46, center = TRUE, scale. = TRUE)
pca_DAM_df_BA46 <- as.data.frame(pca_DAM_BA46$x)


rosnerTest(pca_DAM_df_BA46$PC1, k = 5)
# no outliers, continue

pca_DAM_df_BA46 <- pca_DAM_df_BA46 %>%
  rownames_to_column(var = "library_name") %>%
  dplyr::select(library_name, PC1, PC2, PC3) %>%
  left_join(sample_meta, by = "library_name")


pca_DAM_df_BA46$Classifcation <- factor(pca_DAM_df_BA46$Classifcation, levels = c("Control", "ALS", "ALS - borderline", "ALSci"))

#########################

sdev_BA46 <- pca_DAM_BA46$sdev
variance_BA46 <- sdev_BA46^2
variance_explained_BA46 <- variance_BA46 / sum(variance_BA46) * 100
# PC1: 21.62%

#########################

# save
saveRDS(pca_DAM_df_BA46, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_DAM_final_BA46.rds")
saveRDS(variance_explained_BA46, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_pct_variance_explained_BA46.rds")

##################################################
##################################################
##################################################

# run GLMs

pca_DAM_df_BA44 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_DAM_final_BA44.rds")
pca_DAM_df_BA46 <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_DAM_final_BA46.rds")

pca_DAM_df_BA44$age_scaled <- scale(pca_DAM_df_BA44$Age)
pca_DAM_df_BA46$age_scaled <- scale(pca_DAM_df_BA46$Age)

##################################################

glm_BA44 <- glm(formula = PC1 ~ Classifcation + Sex + age_scaled, data = pca_DAM_df_BA44)
model_summary_BA44 <- summary(glm_BA44)
coef_BA44 <- model_summary_BA44$coefficients
coef_BA44 <- coef_BA44 %>%
  as.data.frame(.) %>%
  rownames_to_column(var = "test")
coef_BA44$test <- paste0(coef_BA44$test, "_BA44")

  
glm_BA46 <- glm(formula = PC1 ~ Classifcation + Sex + age_scaled, data = pca_DAM_df_BA46)
model_summary_BA46 <- summary(glm_BA46)
coef_BA46 <- model_summary_BA46$coefficients
coef_BA46 <- coef_BA46 %>%
  as.data.frame(.) %>%
  rownames_to_column(var = "test")
coef_BA46$test <- paste0(coef_BA46$test, "_BA44")

#########################

# merge

glm_summary_merged <- as.data.frame(rbind(coef_BA44, coef_BA46))


glm_summary_merged <- glm_summary_merged[grepl("Class", glm_summary_merged$test),]
glm_summary_merged$dataset <- "Petrescu_2025"

rownames(glm_summary_merged) <- NULL

#########################

# save
saveRDS(glm_summary_merged, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Petrescu_2025/PCA_GLM_summary.rds")

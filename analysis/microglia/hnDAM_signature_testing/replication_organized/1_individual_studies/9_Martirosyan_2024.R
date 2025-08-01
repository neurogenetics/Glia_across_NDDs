library(edgeR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(SeuratDisk)
library(EnvStats)
library(Matrix)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

cts <- read.csv("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/processed_data_from_study/GSE243639_Filtered_count_table.csv")
cts <- cts %>%
  column_to_rownames(var = "X")
cts <- as(as.matrix(cts), "dgCMatrix")

seurat <- CreateSeuratObject(counts = cts)

##################################################

# dataset is already filtered, don't need to QC just process

# add donor to metadata
seurat$donor <- sub("_.*", "", colnames(seurat))


seurat$pct.mito <- PercentageFeatureSet(seurat, pattern = "^MT")


seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat, vars.to.regress = "pct.mito")
seurat <- RunPCA(seurat)

seurat <- RunHarmony(object = seurat, reduction = "pca", group.by.vars = "donor", 
                              reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "RNA")

ElbowPlot(seurat, ndims = 50, reduction = "harmony")
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:15)
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:15)
seurat <- FindClusters(seurat, resolution = 0.1)


DimPlot_scCustom(seurat)

DotPlot_scCustom(seurat, features = c("C3", "CSF1R", "RBFOX3", "AQP4", "PLP1", "PDGFRA"))

# microglia are cluster 2

micro <- subset(seurat, cells = colnames(seurat)[seurat$seurat_clusters == "2"])

saveRDS(micro, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_subset.rds")

##################################################

Martirosyan_micro <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_subset.rds")


# pseudobulk by donor
pb_cts <- as.data.frame(AggregateExpression(Martirosyan_micro, assays = "RNA", group.by = "donor"))
colnames(pb_cts) <- gsub("RNA.", "", colnames(pb_cts))


# retain samples w/ at least 50 cells
t <- Martirosyan_micro@meta.data %>%
  group_by(donor) %>%
  summarise(n = n())

donors_50cells <- t$donor[t$n >= 50]

pb_cts <- pb_cts[,colnames(pb_cts) %in% donors_50cells]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(Martirosyan_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_pb_cpm.rds")
saveRDS(Martirosyan_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

donor_meta <- read.csv("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/Martirosyan_2024_donor_meta.csv")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_pb_cpm.rds")

##################################################

donor_meta <- donor_meta[donor_meta$Sample.ID %in% colnames(pb_cpm), ]

donor_meta <- donor_meta %>%
  dplyr::select(Sample.ID, Clinical.diagnosis, Age, Sex, PMI.hours, RIN.measure) %>%
  dplyr::rename(donor = Sample.ID, group = Clinical.diagnosis, age = Age, sex = Sex, pmi = PMI.hours, rin = RIN.measure)


donor_meta$pmi_scaled <- scale(donor_meta$pmi)
donor_meta$age_scaled <- scale(donor_meta$age)

donor_meta$group <- factor(donor_meta$group, levels = c("HC", "PD"))


# save
saveRDS(donor_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/donor_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
donor_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/donor_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])
rownames(pb_cpm_DAM) <- gsub("_.*", "", rownames(pb_cpm_DAM))


# filter for usable donors
pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donor_meta$donor,]

pb_cpm_DAM <- pb_cpm_DAM[, !colnames(pb_cpm_DAM) == "CTSD"]

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


# 1 outlier, remove & re-run PCA
donors_keep <- pca_DAM_df$donor[pca_DAM_df$PC1 > -18]

pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donors_keep,]

pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
pca_DAM_df <- as.data.frame(pca_DAM$x)

pca_DAM_df <- pca_DAM_df %>%
  rownames_to_column(var = "donor") %>%
  dplyr::select(donor, PC1, PC1, PC2, PC3) %>%
  left_join(donor_meta, by = "donor")

rosnerTest(pca_DAM_df$PC1, k = 10)


# 1 outlier, remove & re-run PCA
donors_keep <- pca_DAM_df$donor[pca_DAM_df$PC1 > -17]

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
# PC1: 26.15%

#########################

# save
saveRDS(pca_DAM_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/PCA_DAM_final.rds")
saveRDS(variance_explained, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/PCA_pct_variance_explained.rds")

##################################################
##################################################
##################################################

# run GLM

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/PCA_DAM_final.rds")

##################################################

glm <- glm(formula = PC1 ~ group + sex + age_scaled + pmi_scaled, data = pca_DAM_df)
model_summary <- summary(glm)
coef <- model_summary$coefficients

#########################

# merge

glm_summary <- coef %>%
  as.data.frame(.) %>%
  rownames_to_column(var = "test")
glm_summary <- glm_summary[grepl("group", glm_summary$test),]
glm_summary$dataset <- "Martirosyan_2024"

rownames(glm_summary) <- NULL

#########################

# save
saveRDS(glm_summary, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/PCA_GLM_summary.rds")

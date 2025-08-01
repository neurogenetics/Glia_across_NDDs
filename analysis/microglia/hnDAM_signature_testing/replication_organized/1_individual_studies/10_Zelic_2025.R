library(edgeR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(SeuratDisk)
library(EnvStats)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# Processed data from Zelic et al., 2025, Immunity, downloaded from GEO (GSE287257)

# read in / merge sample objects

samples <- list.files("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/processed_data_from_study/", recursive = F, full.names = F)
# sample = samples[1]

seurat_list <- list()

for (sample in samples){
  sample <- gsub("_filtered_feature_bc_matrix.h5", "", sample)
  print(paste0("Processing ", sample))
  cts <- Read10X_h5(paste0("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/processed_data_from_study/", sample, "_filtered_feature_bc_matrix.h5"))
  seurat <- CreateSeuratObject(counts = cts)
  
  seurat$sample <- sample
  colnames(seurat) <- paste0(colnames(seurat), "_", sample)
  
  seurat_list[[sample]] <- seurat
}

seurat_merged <- Merge_Seurat_List(seurat_list)

seurat_merged[["RNA"]] <- JoinLayers(seurat_merged[["RNA"]])

##################################################

# process data

# QC

seurat_merged$pct.mito <- PercentageFeatureSet(seurat_merged, pattern = "^MT")

# same filters as used in the paper
seurat_filtered <- subset(seurat_merged, nFeature_RNA > 1000 & nCount_RNA > 2500 & pct.mito < 15)


seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- FindVariableFeatures(seurat_filtered)
seurat_filtered <- ScaleData(seurat_filtered, vars.to.regress = "pct.mito")
seurat_filtered <- RunPCA(seurat_filtered)

seurat_filtered <- RunHarmony(object = seurat_filtered, reduction = "pca", group.by.vars = "sample", 
                              reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "RNA")

ElbowPlot(seurat_filtered, ndims = 50, reduction = "harmony")
seurat_filtered <- FindNeighbors(seurat_filtered, reduction = "harmony", dims = 1:30)
seurat_filtered <- RunUMAP(seurat_filtered, reduction = "harmony", dims = 1:30)
seurat_filtered <- FindClusters(seurat_filtered, resolution = 0.1)


DimPlot_scCustom(seurat_filtered)

DotPlot_scCustom(seurat_filtered, features = c("C3", "CSF1R", "RBFOX3", "AQP4", "PLP1", "PDGFRA"))

# microglia are cluster 1

micro <- subset(seurat_filtered, cells = colnames(seurat_filtered)[seurat_filtered$seurat_clusters == "1"])

saveRDS(micro, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_subset_unprocessed.rds")

##################################################

Zelic_micro <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_subset_unprocessed.rds")


# pseudobulk by donor
pb_cts <- as.data.frame(AggregateExpression(Zelic_micro, assays = "RNA", group.by = "sample"))
colnames(pb_cts) <- gsub("RNA.", "", colnames(pb_cts))
colnames(pb_cts) <- gsub("\\.", "_", colnames(pb_cts))


# retain samples w/ at least 50 cells
t <- Zelic_micro@meta.data %>%
  group_by(sample) %>%
  summarise(n = n())

samples_50cells <- t$sample[t$n >= 50]

pb_cts <- pb_cts[,colnames(pb_cts) %in% samples_50cells]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(Zelic_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_pb_cpm.rds")
saveRDS(Zelic_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

donor_meta <- read.csv("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/Zelic_2025_donor_meta.csv")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_pb_cpm.rds")

##################################################

donor_meta$group <- factor(donor_meta$group, levels = c("HC", "ALS"))

donor_meta$age_scaled <- scale(donor_meta$age)
donor_meta$pmi_scaled <- scale(donor_meta$pmi)


# save
saveRDS(donor_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/donor_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
donor_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/donor_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])
rownames(pb_cpm_DAM) <- gsub("^[^_]*_", "", rownames(pb_cpm_DAM))


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
# no more outliers, proceed

#########################

sdev <- pca_DAM$sdev
variance <- sdev^2
variance_explained <- variance / sum(variance) * 100
# PC1: 39.08%

#########################

# save
saveRDS(pca_DAM_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/PCA_DAM_final.rds")
saveRDS(variance_explained, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/PCA_pct_variance_explained.rds")

##################################################
##################################################
##################################################

# run GLM

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/PCA_DAM_final.rds")

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
glm_summary$dataset <- "Zelic_2025"

rownames(glm_summary) <- NULL

#########################

# save
saveRDS(glm_summary, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/PCA_GLM_summary.rds")

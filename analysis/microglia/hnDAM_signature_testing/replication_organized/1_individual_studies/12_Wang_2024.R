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

# Processed data from Wang et al., 2024, Sci. Adv.
# downloaded from Human Cell Atlas (https://explore.data.humancellatlas.org/projects/16cd6791-2adb-4d0f-8222-0184dada6456)

# read in / merge sample objects

samples <- list.dirs("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/processed_data_from_study/", recursive = F, full.names = F)
# sample = samples[1]

seurat_list <- list()

for (sample in samples){
  print(paste0("Processing ", sample))
  cts <- Read10X(paste0("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/processed_data_from_study/", sample, "/filtered_feature_bc_matrix/"))
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

# VlnPlot_scCustom(seurat_merged, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0)

seurat_filtered <- subset(seurat_merged, nFeature_RNA < 2500 & nFeature_RNA > 200)

seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- FindVariableFeatures(seurat_filtered)
seurat_filtered <- ScaleData(seurat_filtered, vars.to.regress = "pct.mito")
seurat_filtered <- RunPCA(seurat_filtered)

seurat_filtered <- RunHarmony(object = seurat_filtered, reduction = "pca", group.by.vars = "sample", 
                              reduction.save = 'harmony', plot_convergence = T, lambda = NULL, assay.use = "RNA")

ElbowPlot(seurat_filtered, ndims = 50, reduction = "harmony")
seurat_filtered <- FindNeighbors(seurat_filtered, reduction = "harmony", dims = 1:15)
seurat_filtered <- RunUMAP(seurat_filtered, reduction = "harmony", dims = 1:15)
seurat_filtered <- FindClusters(seurat_filtered, resolution = 0.1)


DimPlot_scCustom(seurat_filtered)

DotPlot_scCustom(seurat_filtered, features = c("C3", "CSF1R", "RBFOX3", "AQP4", "PLP1", "PDGFRA"))

# microglia are cluster 1

micro <- subset(seurat_filtered, cells = colnames(seurat_filtered)[seurat_filtered$seurat_clusters == "1"])

saveRDS(micro, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_subset_unprocessed.rds")

##################################################

Wang_micro <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_subset_unprocessed.rds")


# pseudobulk by donor
pb_cts <- as.data.frame(AggregateExpression(Wang_micro, assays = "RNA", group.by = "sample"))
colnames(pb_cts) <- gsub("RNA.", "", colnames(pb_cts))
colnames(pb_cts) <- gsub("\\.", "_", colnames(pb_cts))


# retain samples w/ at least 50 cells
t <- Wang_micro@meta.data %>%
  group_by(sample) %>%
  summarise(n = n())

samples_50cells <- t$sample[t$n >= 50]

pb_cts <- pb_cts[,colnames(pb_cts) %in% samples_50cells]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(Wang_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_pb_cpm.rds")
saveRDS(Wang_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

donor_meta <- read.csv("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/wang_2024_donor_meta.csv")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_pb_cpm.rds")

##################################################

colnames(pb_cpm) <- gsub("_.*", "", colnames(pb_cpm))

donor_meta <- donor_meta[donor_meta$donor %in% colnames(pb_cpm), ]

donor_meta$pmi_scaled <- scale(donor_meta$pmi)
donor_meta$age_scaled <- scale(donor_meta$age_at_death)

donor_meta$group <- ifelse(donor_meta$group %in% c("PD", "PDD"), "PD", "HC")
donor_meta$group <- factor(donor_meta$group, levels = c("HC", "PD"))


# save
saveRDS(donor_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/donor_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
donor_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/donor_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])
rownames(pb_cpm_DAM) <- gsub("_.*", "", rownames(pb_cpm_DAM))


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


# 1 outlier, remove & re-run PCA
donors_keep <- pca_DAM_df$donor[pca_DAM_df$PC1 < 16]

pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donors_keep,]

pb_cpm_DAM <- pb_cpm_DAM[, !colnames(pb_cpm_DAM) == "CTSD"]

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
# PC1: 25.02%

#########################

# save
saveRDS(pca_DAM_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/PCA_DAM_final.rds")
saveRDS(variance_explained, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/PCA_pct_variance_explained.rds")

##################################################
##################################################
##################################################

# run GLM

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/PCA_DAM_final.rds")

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
glm_summary$dataset <- "Wang_2024"

rownames(glm_summary) <- NULL

#########################

# save
saveRDS(glm_summary, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Wang_2024/PCA_GLM_summary.rds")

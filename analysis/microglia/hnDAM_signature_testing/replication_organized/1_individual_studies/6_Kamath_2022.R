library(edgeR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(SeuratDisk)
library(EnvStats)
library(harmony)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################

# prep data & counts for PCA analysis

##################################################

# processed data downloaded from Broad institute single-cell portal (https://singlecell.broadinstitute.org/single_cell/study/SCP1768/)
cts <- Read10X(data.dir = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/processed_data_from_study/")


# this is the file called micro_UMAP 
micro_meta <- read.csv("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/Kamath_microglia.csv")


meta <- read.csv("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/processed_data_from_study/metadata.tsv", sep = "\t")
meta <- meta[-1, ]
rownames(meta) <- NULL
meta <- meta %>%
  column_to_rownames(var = "NAME")


meta <- meta[rownames(meta) %in% micro_meta$NAME, ]
meta <- meta[meta$FACS_Classification == "Negative", ]

micro_cts <- cts[, colnames(cts) %in% rownames(meta)]

Kamath_micro <- CreateSeuratObject(micro_cts, meta.data = meta)

##################################################

# pseudobulk by donor
pb_cts <- as.data.frame(AggregateExpression(Kamath_micro, assays = "RNA", group.by = "donor_id"))
colnames(pb_cts) <- gsub("RNA.g", "", colnames(pb_cts))


# retain donors w/ at least 50 cells
t <- Kamath_micro@meta.data %>%
  group_by(donor_id) %>%
  summarise(n = n())

donors_50cells <- t$donor_id[t$n >= 50]

pb_cts <- pb_cts[,colnames(pb_cts) %in% donors_50cells]


# TMM-normalize counts
pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))


# get single-cell counts for later
sc_cts <- LayerData(Kamath_micro, assay = "RNA", layer = "counts")


# save
saveRDS(pb_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/microglia_pb_cts.rds")
saveRDS(pb_cpm, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/microglia_pb_cpm.rds")
saveRDS(Kamath_micro@meta.data, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/microglia_celllevel_meta.rds")
saveRDS(sc_cts, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/microglia_sc_cts.rds")

##################################################
##################################################
##################################################

# prep metadata and filter for donors to keep 

celllevel_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/microglia_celllevel_meta.rds")
pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/microglia_pb_cpm.rds")

##################################################

donor_meta <- celllevel_meta %>%
  distinct(donor_id, .keep_all = T)

rownames(donor_meta) <- NULL

donor_meta <- donor_meta %>%
  dplyr::select(donor_id, disease__ontology_label, Donor_Age, Donor_PMI, sex) %>%
  dplyr::rename(donor = donor_id, group = disease__ontology_label, age = Donor_Age, pmi = Donor_PMI)

donor_meta$group <- ifelse(donor_meta$group == "normal", "HC", "PD/DLB")

donor_meta$group <- factor(donor_meta$group, levels = c("HC", "PD/DLB"))


donor_meta$pmi <- as.numeric(donor_meta$pmi)
donor_meta$age <- as.numeric(donor_meta$age)


donor_meta$pmi_scaled <- scale(donor_meta$pmi)
donor_meta$age_scaled <- scale(donor_meta$age)


# save
saveRDS(donor_meta, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/donor_meta.rds")

##################################################
##################################################
##################################################

# run PCA analysis

pb_cpm <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/microglia_pb_cpm.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
donor_meta <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/donor_meta.rds")

##################################################

# filter for DAM signature genes
pb_cpm_DAM <- t(pb_cpm[rownames(pb_cpm) %in% c(DAM_genes), ])
rownames(pb_cpm_DAM) <- gsub("_.*", "", rownames(pb_cpm_DAM))


# filter for usable donors
pb_cpm_DAM <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donor_meta$donor,]

pb_cpm_DAM <- pb_cpm_DAM[, !colnames(pb_cpm_DAM) %in% c("RAB7B", "LINC01235")]


# run PCA
pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
pca_DAM_df <- as.data.frame(pca_DAM$x)


# join PCA df w/ metadata
pca_DAM_df <- pca_DAM_df %>%
  rownames_to_column(var = "donor") %>%
  dplyr::select(donor, PC1, PC2, PC3) %>%
  left_join(donor_meta, by = "donor")


# test for outliers
rosnerTest(pca_DAM_df$PC1, k = 5)
# no outliers, proceed

#########################

sdev <- pca_DAM$sdev
variance <- sdev^2
variance_explained <- variance / sum(variance) * 100
# PC1: 22.09%

#########################

# save
saveRDS(pca_DAM_df, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/PCA_DAM_final.rds")
saveRDS(variance_explained, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/PCA_pct_variance_explained.rds")

##################################################
##################################################
##################################################

# run GLM

pca_DAM_df <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/PCA_DAM_final.rds")

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
glm_summary$dataset <- "Kamath_2022"

rownames(glm_summary) <- NULL

#########################

# save
saveRDS(glm_summary, file = "./analysis/microglia/DAM_signature_replication/individual_studies/Kamath_2022/PCA_GLM_summary.rds")

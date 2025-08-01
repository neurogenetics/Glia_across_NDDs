library(edgeR)
library(tidyverse)
library(ggplot2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

pb_cts <- readRDS("./analysis/microglia/differential_expression/cts_pseudobulked_donor_x_region.rds")

DAM_genes_up <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")

diffexp_tests <- read.delim("./analysis/diffexp_tests.txt", header = F)
diffexp_tests <- diffexp_tests$V1

donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")
donor_meta$group <- gsub("-", "", donor_meta$group)

##################################################
##################################################
##################################################

# TMM-normalize counts

pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))

##################################################

# test = diffexp_tests[1]

PCA_DAM_df_list <- list()
GLM_res_summary <- data.frame()
variance_explained_list <- list()

for (test in diffexp_tests){
  components <- unlist(strsplit(test, "_"))
  dataset <- components[1]
  brain_region <- components[2]
  disease <- components[3]
  rm(components)

  # filter for donors in this comparison
  donor_meta_filtered <- donor_meta[donor_meta$study == dataset & donor_meta$group %in% c("HC", disease), ]
  donor_meta_filtered$donor_region <- paste0(donor_meta_filtered$donor, "_", brain_region)
  
  # filter cts matrix for donors and DAM genes
  pb_cpm_filtered <- pb_cpm[, colnames(pb_cpm) %in% donor_meta_filtered$donor_region]
  pb_cpm_DAM <- t(pb_cpm_filtered[rownames(pb_cpm_filtered) %in% DAM_genes_up, ])
  pb_cpm_DAM <- pb_cpm_DAM[, colSums(pb_cpm_DAM != 0) > 0]
  
  # run PCA
  pca_DAM <- prcomp(pb_cpm_DAM, center = TRUE, scale. = TRUE)
  
  # organize PCA results w/ donor metadata
  pca_DAM_df <- as.data.frame(pca_DAM$x)
  rownames(pca_DAM_df) <- gsub(paste0("_", brain_region), "", rownames(pca_DAM_df))
  pca_DAM_df <- pca_DAM_df %>%
    rownames_to_column(var = "donor") %>%
    dplyr::select(donor, PC1, PC2, PC3) %>%
    left_join(donor_meta_filtered, by = "donor")
  pca_DAM_df$age_scaled <- scale(pca_DAM_df$age)
  pca_DAM_df$group <- factor(pca_DAM_df$group, levels = c("HC", disease))
  
  # get % variance explained by PCs
  sdev <- pca_DAM$sdev
  variance <- sdev^2
  variance_explained <- variance / sum(variance) * 100
  
  # run the GLM
  glm <- glm(formula = PC1 ~ group + sex + age_scaled, data = pca_DAM_df)
  model_summary <- summary(glm)
  coef <- model_summary$coefficients
  coef <- coef %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "test")
  coef$test <- paste0(coef$test, "_", brain_region)
  
  GLM_res_summary <- rbind(GLM_res_summary, coef)
  PCA_DAM_df_list[[test]] <- pca_DAM_df
  variance_explained_list[[test]] <- variance_explained
}

GLM_res_summary <- GLM_res_summary[grepl("group", GLM_res_summary$test), ]
GLM_res_summary$padj <- p.adjust(GLM_res_summary$`Pr(>|t|)`, method = "BH")

#########################

saveRDS(GLM_res_summary, file = "./analysis/microglia/DAM_signature_discovery/PCA_GLM_summary.rds")
saveRDS(PCA_DAM_df_list, file = "./analysis/microglia/DAM_signature_discovery/PCA_df_list.rds")
saveRDS(variance_explained_list, file = "./analysis/microglia/DAM_signature_discovery/variance_explained_list.rds")

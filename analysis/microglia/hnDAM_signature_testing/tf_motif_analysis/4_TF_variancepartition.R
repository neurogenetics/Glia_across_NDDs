library(edgeR)
library(tidyverse)
library(ggplot2)
library(vegan)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# variance partitioning analysis w/ vegan

pb_cts <- readRDS("./analysis/microglia/differential_expression/cts_pseudobulked_donor_x_region.rds")

DAM_genes_up <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")

donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")
donor_meta$group <- gsub("-", "", donor_meta$group)

##################################################

# TMM-normalize counts

pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(t(edgeR::cpm(pb_dge)))

##################################################

tf_expr <- data.frame(sample = rownames(pb_cpm),
                      MITF = pb_cpm$MITF,
                      ARID5B = pb_cpm$ARID5B,
                      MAFB = pb_cpm$MAFB,
                      PPARG = pb_cpm$PPARG,
                      SPI1 = pb_cpm$SPI1,
                      CEBPA = pb_cpm$CEBPA)
tf_expr$donor = sub("_(?!.*_).*", "", tf_expr$sample, perl = TRUE)
tf_expr$region = sub(".*_(.*)", "\\1", tf_expr$sample)

pb_cpm_DAM <- pb_cpm[, colnames(pb_cpm) %in% DAM_genes_up]

##################################################

tests <- c("Mathys_HIP", "Mathys_EC", "Mathys_TH", "Mathys_FC/PFC", "Mathys_AnG", "Mathys_TC/MTG",
           "AMP-PD_DMV", "AMP-PD_GPi", "AMP-PD_M1", "AMP-PD_FC/PFC", "AMP-PD_OC/V1",
           "Gerrits_FC/PFC", "Gerrits_TC/MTG", "Gerrits_OC/V1",
           "Pineda_M1", "Pineda_FC/PFC")

# test = tests[1]

##################################################
##################################################
##################################################

# for the four TFs in the module

summary_unique_all <- data.frame()
summary_contr_all <- data.frame()
anova_res_all <- data.frame()

for (test in tests){
  components <- unlist(strsplit(test, "_"))
  dataset <- components[1]
  brain_region <- components[2]

  # filter for donors in this comparison
  donor_meta_filtered <- donor_meta[donor_meta$study == dataset, ]
  donor_meta_filtered$donor_region <- paste0(donor_meta_filtered$donor, "_", brain_region)
  
  # filter cts matrix for donors and DAM genes
  pb_cpm_filtered <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donor_meta_filtered$donor_region, ]
  pb_cpm_filtered <- pb_cpm_filtered[, colSums(pb_cpm_filtered != 0) > 0]
  pb_cpm_filtered <- pb_cpm_filtered %>%
    dplyr::select(-MITF, -ARID5B, -MAFB, -PPARG)
  
  # run PCA
  pca_DAM <- prcomp(pb_cpm_filtered, center = TRUE, scale. = TRUE)
  pca_DAM_df <- as.data.frame(pca_DAM$x)
  rownames(pca_DAM_df) <- gsub(paste0("_", brain_region), "", rownames(pca_DAM_df))
  pca_DAM_df <- pca_DAM_df %>%
    rownames_to_column(var = "donor") %>%
    dplyr::select(donor, PC1, PC2, PC3) %>%
    left_join(donor_meta_filtered, by = "donor")
  pca_DAM_df <- na.omit(pca_DAM_df)
  pca_DAM_df$age_scaled <- scale(pca_DAM_df$age)

  tf_expr_filtered <- tf_expr[tf_expr$region == brain_region, ]
  pca_DAM_df <- pca_DAM_df %>%
    left_join(tf_expr_filtered, by = "donor")
  pca_DAM_df$MITF <- scale(pca_DAM_df$MITF)
  pca_DAM_df$MAFB <- scale(pca_DAM_df$MAFB)
  pca_DAM_df$ARID5B <- scale(pca_DAM_df$ARID5B)
  pca_DAM_df$PPARG <- scale(pca_DAM_df$PPARG)
  
  
  other_vars <- subset(pca_DAM_df, select = c(sex, age_scaled))
  tf_vars <- subset(pca_DAM_df, select = c(ARID5B, MAFB, MITF, PPARG))
  
  resid_pc1 <- resid(lm(PC1 ~ age_scaled + sex, data = pca_DAM_df))
  
  part_tf_only <- varpart(resid_pc1,
                          pca_DAM_df["ARID5B"],
                          pca_DAM_df["MAFB"],
                          pca_DAM_df["MITF"],
                          pca_DAM_df["PPARG"])
  summary <- summary(part_tf_only)
  
  summary_unique <- as.data.frame(summary$uniqpart)
  summary_contr <- as.data.frame(summary$contribpart)
  
  rownames(summary_unique) <- NULL
  rownames(summary_contr) <- NULL
  
  summary_unique$TF <- c("ARID5B", "MAFB", "MITF", "PPARG")
  summary_contr$TF <- c("ARID5B", "MAFB", "MITF", "PPARG")
  
  summary_unique$dataset_region <- test
  summary_contr$dataset_region <- test
  
  summary_unique_all <- rbind(summary_unique_all, summary_unique)
  summary_contr_all <- rbind(summary_contr_all, summary_contr)
  
  anova_1 <- anova.cca(rda(pca_DAM_df$PC1 ~ ARID5B + Condition(MAFB + MITF + PPARG + age_scaled + sex), data = pca_DAM_df), permutations = 10000)
  anova_2 <- anova.cca(rda(pca_DAM_df$PC1 ~ MAFB + Condition(ARID5B + MITF + PPARG + age_scaled + sex), data = pca_DAM_df), permutations = 10000)
  anova_3 <- anova.cca(rda(pca_DAM_df$PC1 ~ MITF + Condition(ARID5B + MAFB + PPARG + age_scaled + sex), data = pca_DAM_df), permutations = 10000)
  anova_4 <- anova.cca(rda(pca_DAM_df$PC1 ~ PPARG + Condition(ARID5B + MAFB + MITF+ age_scaled + sex), data = pca_DAM_df), permutations = 10000)
  anova_5 <- anova.cca(rda(pca_DAM_df$PC1 ~ ARID5B + MAFB + MITF + PPARG + Condition(age_scaled + sex), data = pca_DAM_df), permutations = 10000)
  
  anova_res <- data.frame(ARID5B = anova_1$`Pr(>F)`[1],
                          MAFB = anova_2$`Pr(>F)`[1], 
                          MITF = anova_3$`Pr(>F)`[1],
                          PPARG = anova_4$`Pr(>F)`[1],
                          combined = anova_5$`Pr(>F)`[1])
  anova_res$dataset_region <- test
  
  anova_res_all <- rbind(anova_res_all, anova_res)
}


anova_res_long <- anova_res_all %>%
  pivot_longer(cols = c("MITF", "ARID5B", "MAFB", "PPARG", "combined"), names_to = "TF") %>%
  dplyr::rename(pval = value)



summary_contr_all[summary_contr_all < 0] <- 0

total_var_explained <- summary_contr_all %>%
  group_by(dataset_region) %>%
  summarise(total_var = sum(`summary$contribpart`))

total_var_explained$TF <- "combined"



colnames(summary_unique_all) <- c("var_explained", "TF", "dataset_region")
colnames(total_var_explained) <- c("dataset_region", "var_explained", "TF")

summary_final <- rbind(summary_unique_all, total_var_explained)


summary_final <- summary_final %>%
  left_join(anova_res_long, by = c("TF", "dataset_region"))

summary_final$pval <- ifelse(summary_final$var_explained < 0, NA, summary_final$pval)

summary_final$var_explained[summary_final$var_explained < 0] <- 0

##################################################

saveRDS(summary_final, file = "./analysis/microglia/DAM_signature_discovery/tf_varpart/final_summary_DAM_TFs.rds")

##################################################
##################################################
##################################################

# for alternate "negative control" TFs

overall_summary <- data.frame()

for (test in tests){
  components <- unlist(strsplit(test, "_"))
  dataset <- components[1]
  brain_region <- components[2]
  
  # filter for donors in this comparison
  donor_meta_filtered <- donor_meta[donor_meta$study == dataset, ]
  donor_meta_filtered$donor_region <- paste0(donor_meta_filtered$donor, "_", brain_region)
  
  # filter cts matrix for donors and DAM genes
  pb_cpm_filtered <- pb_cpm_DAM[rownames(pb_cpm_DAM) %in% donor_meta_filtered$donor_region, ]
  pb_cpm_filtered <- pb_cpm_filtered[, colSums(pb_cpm_filtered != 0) > 0]
  pb_cpm_filtered <- pb_cpm_filtered %>%
    dplyr::select(-MITF, -ARID5B, -MAFB, -PPARG)
  
  # run PCA
  pca_DAM <- prcomp(pb_cpm_filtered, center = TRUE, scale. = TRUE)
  pca_DAM_df <- as.data.frame(pca_DAM$x)
  rownames(pca_DAM_df) <- gsub(paste0("_", brain_region), "", rownames(pca_DAM_df))
  pca_DAM_df <- pca_DAM_df %>%
    rownames_to_column(var = "donor") %>%
    dplyr::select(donor, PC1, PC2, PC3) %>%
    left_join(donor_meta_filtered, by = "donor")
  pca_DAM_df <- na.omit(pca_DAM_df)
  pca_DAM_df$age_scaled <- scale(pca_DAM_df$age)
  
  tf_expr_filtered <- tf_expr[tf_expr$region == brain_region, ]
  pca_DAM_df <- pca_DAM_df %>%
    left_join(tf_expr_filtered, by = "donor")
  pca_DAM_df$MITF <- scale(pca_DAM_df$MITF)
  pca_DAM_df$MAFB <- scale(pca_DAM_df$MAFB)
  pca_DAM_df$ARID5B <- scale(pca_DAM_df$ARID5B)
  pca_DAM_df$PPARG <- scale(pca_DAM_df$PPARG)
  
  # residualize on age/sex
  resid_pc1 <- resid(lm(PC1 ~ age_scaled + sex, data = pca_DAM_df))
  DAM_tf_vars <- subset(pca_DAM_df, select = c(MITF, MAFB, ARID5B, PPARG))
  
  # loop thru alt TFs to test, run the varpart using the alt TF expression + a df w/ the DAM tf expressions
  # condition the anova on the expression of the DAM TFs
  for (tf in c("SPI1", "CEBPA")){
    pca_DAM_df[[tf]] <- scale(pca_DAM_df[[tf]])
    
    varpart <- varpart(resid_pc1, pca_DAM_df[[tf]], DAM_tf_vars)
    summary <- summary(varpart)
    summary_unique <- as.data.frame(summary$uniqpart)

    anova_1 <- anova.cca(rda(pca_DAM_df$PC1 ~ pca_DAM_df[[tf]] + Condition(MITF + ARID5B + PPARG + MAFB + age_scaled + sex), data = pca_DAM_df), permutations = 10000)
    
    summary <- data.frame(tf_var_expl = summary_unique[1,1], 
                          DAM_tf_var_expl = summary_unique[2,1],
                          TF = tf,
                          anova = anova_1[1,4], 
                          test = test)
    
    overall_summary <- rbind(overall_summary, summary)
    }
}


overall_summary$anova[overall_summary$tf_var_expl < 0] <- NA

overall_summary[overall_summary < 0] <- 0
overall_summary[overall_summary < 0] <- 0

##################################################

saveRDS(overall_summary, file = "./analysis/microglia/DAM_signature_discovery/tf_varpart/other_TF_contributions_summary.rds")

library(tidyverse)
library(mashr)
library(DESeq2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

source("./code_organized/functions/rename_datasets.R")

##################################################
##################################################
##################################################

# microglia markers

t <- readRDS("./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea_ANNOTATED.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/all_microglia_markers.csv", row.names = F)

##################################################

# ALL astrocyte markers

t <- readRDS("./analysis/astrocytes/cluster_characterization/all_astros_markers_ALL_ANNOTATED.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/all_astrocyte_markers.csv", row.names = F)

##################################################

# cortical astrocyte markers

t <- readRDS("./analysis/astrocytes/cluster_characterization/cortical_astros_markers_ALL_ANNOTATED.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/all_cortical_astrocyte_markers.csv", row.names = F)

##################################################

# cortical astrocyte markers

t <- readRDS("./analysis/astrocytes/cluster_characterization/subcortical_astros_markers_ALL_ANNOTATED.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/all_subcortical_astrocyte_markers.csv", row.names = F)

##################################################

# oligodendrocyte markers

t <- readRDS("./analysis/oligodendrocytes/cluster_characterization/oligos_markers_ALL_ANNOTATED.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/all_oligos_markers.csv", row.names = F)

##################################################

# cNMF GEPs, Gerrits

t <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_gene_scores.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/Gerrits_GEPs_gene_scores.csv", row.names = T)

##################################################

# cNMF GEPs, Mathys

t <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_gene_scores.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/Mathys_GEPs_gene_scores.csv", row.names = T)

##################################################

# cNMF GEPs, Pineda

t <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_gene_scores.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/Pineda_GEPs_gene_scores.csv", row.names = T)

##################################################

# cNMF GEPs, AMPPD

t <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_gene_scores.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/AMPPD_GEPs_gene_scores.csv", row.names = T)

##################################################

# DAM sig GOBP

t <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/DAM_sig_GOBP_all.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/DAM_sig_GOBP_all.csv", row.names = T)

##################################################

# proportions summary stats (from GLMs)

t1 <- readRDS("./analysis/microglia/cluster_proportions/micro_props_logit_GLM_results_1v1.rds")
t2 <- readRDS("./analysis/astrocytes/cluster_proportions/cortical_astro_props_logit_GLM_results_1v1.rds")
t3 <- readRDS("./analysis/oligodendrocytes/cluster_proportions/oligo_props_logit_GLM_results_1v1.rds")

glmres <- rbind(t1, t2, t3)
colnames(glmres) <- c("cluster", "term", "beta", "stderr", "t-value", "p-value", "dataset_region", "FDR")
glmres <- glmres[c("cluster", "dataset_region", "term", "beta", "stderr", "t-value", "p-value",  "FDR")]

write.csv(glmres, file = "./analysis/csvs_for_supp_tables/proportions_raw_GLM_summarystats.csv", row.names = F)

##################################################

# mashr- corrected GLM summary stats (beta/lfsr)

t1 <- readRDS("./analysis/microglia/cluster_proportions/micro_props_mashr_obj_1v1.rds")
t2 <- readRDS("./analysis/astrocytes/cluster_proportions/cortical_astros_props_mashr_obj_canonical.rds")
t3 <- readRDS("./analysis/oligodendrocytes/cluster_proportions/oligo_props_mashr_obj_canonical.rds")


micro_beta <- as.data.frame(t1$result$PosteriorMean)
micro_lfsr <- as.data.frame(t1$result$lfsr)

micro_lfsr <- micro_lfsr %>%
  tibble::rownames_to_column(var = "cluster") %>%
  pivot_longer(
    cols = -cluster,              
    names_to = "disease_region",  
    values_to = "lfsr") 

micro_beta <- micro_beta %>%
  tibble::rownames_to_column(var = "cluster") %>%
  pivot_longer(
    cols = -cluster,              
    names_to = "disease_region",  
    values_to = "beta_mashr") %>%
  left_join(micro_lfsr, by = c("disease_region", "cluster"))


astro_beta <- as.data.frame(t2$result$PosteriorMean)
astro_lfsr <- as.data.frame(t2$result$lfsr)

astro_lfsr <- astro_lfsr %>%
  tibble::rownames_to_column(var = "cluster") %>%
  pivot_longer(
    cols = -cluster,              
    names_to = "disease_region",  
    values_to = "lfsr") 

astro_beta <- astro_beta %>%
  tibble::rownames_to_column(var = "cluster") %>%
  pivot_longer(
    cols = -cluster,              
    names_to = "disease_region",  
    values_to = "beta_mashr") %>%
  left_join(astro_lfsr, by = c("disease_region", "cluster"))


oligo_beta <- as.data.frame(t3$result$PosteriorMean)
oligo_lfsr <- as.data.frame(t3$result$lfsr)

oligo_lfsr <- oligo_lfsr %>%
  tibble::rownames_to_column(var = "cluster") %>%
  pivot_longer(
    cols = -cluster,              
    names_to = "disease_region",  
    values_to = "lfsr") 

oligo_beta <- oligo_beta %>%
  tibble::rownames_to_column(var = "cluster") %>%
  pivot_longer(
    cols = -cluster,              
    names_to = "disease_region",  
    values_to = "beta_mashr") %>%
  left_join(oligo_lfsr, by = c("disease_region", "cluster"))


mashr_df <- rbind(micro_beta, astro_beta, oligo_beta)

mashr_df$disease_region <- rename_datasets(mashr_df$disease_region)

write.csv(mashr_df, file = "./analysis/csvs_for_supp_tables/proportions_mashr_summarystats.csv", row.names = F)

##################################################

# microglia diffexp

t <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_results_for_mashr.rds")

for (test in names(t)){
  df <- t[[test]]
  colnames(df) <- c("gene", 
                    "nebula_logFC", "nebula_SE", "nebula_pval",
                    "DESeq_basemean", "DESeq_logFC", "DESeq_lfcSE", "DESeq_stat", "DESeq_pval", "DESeq_padj",
                    "nebula_zscore", "DESeq_zscore")
  
  name <- rename_datasets(test)
  
  write.csv(df, file = paste0("./analysis/csvs_for_supp_tables/microglia_diffexp/", name, ".csv"), row.names = F)
}

##################################################

# astrocytes diffexp

t <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_diffexp_results_for_mashr.rds")

for (test in names(t)){
  df <- t[[test]]
  colnames(df) <- c("gene", 
                    "nebula_logFC", "nebula_SE", "nebula_pval",
                    "DESeq_basemean", "DESeq_logFC", "DESeq_lfcSE", "DESeq_stat", "DESeq_pval", "DESeq_padj",
                    "nebula_zscore", "DESeq_zscore")
  
  name <- rename_datasets(test)
  
  write.csv(df, file = paste0("./analysis/csvs_for_supp_tables/astrocytes_diffexp/", name, ".csv"), row.names = F)
}

##################################################

# oligodendrocytes diffexp

t <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_results_for_mashr.rds")

for (test in names(t)){
  df <- t[[test]]
  colnames(df) <- c("gene", 
                    "nebula_logFC", "nebula_SE", "nebula_pval",
                    "DESeq_basemean", "DESeq_logFC", "DESeq_lfcSE", "DESeq_stat", "DESeq_pval", "DESeq_padj",
                    "nebula_zscore", "DESeq_zscore")
  
  name <- rename_datasets(test)
  
  write.csv(df, file = paste0("./analysis/csvs_for_supp_tables/oligodendrocytes_diffexp/", name, ".csv"), row.names = F)
}

##################################################

# DAM replication datasets summary

t1 <- readRDS("./analysis/microglia/DAM_signature_replication/GLM_summary_merged_FDRcorr.rds")
t2 <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_DAM_GSEA_summary.rds")

t1 <- t1 %>%
  left_join(t2, by = "test")

write.csv(t1, file = "./analysis/csvs_for_supp_tables/DAM_replication_summary.csv", row.names = F)

##################################################

# k-means clusters assignments

# microglia

t <- readRDS("./analysis/microglia/differential_expression/microglia_tPCA_kmeans_clusters.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/microglia_kmeans_assignments.csv", row.names = F)


# astros
t <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_tPCA_kmeans_clusters.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/astros_kmeans_assignments.csv", row.names = F)


# oligos
t <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_tPCA_kmeans_clusters.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/oligos_kmeans_assignments.csv", row.names = F)

##################################################

# RUVseq components (only the ones we used)

# microglia

files <- list.files("./analysis/microglia/differential_expression/covar_selection/")

for (file in files){
  covars = readRDS(paste0("./analysis/microglia/differential_expression/covar_selection/", file))
  nRUV = covars$optimal_RUV_W + 1

  RUV_df <- covars$RUV_df
  RUV_df <- RUV_df[,1:nRUV]
  
  name <- gsub("_covar_selection.rds", "", file)
  name <- rename_datasets(name)
  
  write.csv(RUV_df, file = paste0("./analysis/csvs_for_supp_tables/microglia_RUVseq_components/", name, ".csv"), row.names = F)
}


# astrocytes

files <- list.files("./analysis/astrocytes/differential_expression/covar_selection/")

for (file in files){
  covars = readRDS(paste0("./analysis/astrocytes/differential_expression/covar_selection/", file))
  nRUV = covars$optimal_RUV_W + 1
  
  RUV_df <- covars$RUV_df
  RUV_df <- RUV_df[,1:nRUV]
  
  name <- gsub("_covar_selection.rds", "", file)
  name <- rename_datasets(name)
  
  write.csv(RUV_df, file = paste0("./analysis/csvs_for_supp_tables/astrocytes_RUVseq_components/", name, ".csv"), row.names = F)
}


# oligos

files <- list.files("./analysis/oligodendrocytes/differential_expression/covar_selection/")

for (file in files){
  covars = readRDS(paste0("./analysis/oligodendrocytes/differential_expression/covar_selection/", file))
  nRUV = covars$optimal_RUV_W + 1
  
  RUV_df <- covars$RUV_df
  RUV_df <- RUV_df[,1:nRUV]
  
  name <- gsub("_covar_selection.rds", "", file)
  name <- rename_datasets(name)
  
  write.csv(RUV_df, file = paste0("./analysis/csvs_for_supp_tables/oligodendrocytes_RUVseq_components/", name, ".csv"), row.names = F)
}

##################################################

# hnDAM varpart summary stats

t <- as.data.frame(readRDS("./analysis/microglia/DAM_signature_discovery/tf_varpart/anova_res_ALLcomps_df.rds"))
t2 <- as.data.frame(readRDS("./analysis/microglia/DAM_signature_discovery/tf_varpart/variance_explained_ALLcomps_df.rds"))


t <- t %>%
  rownames_to_column(var = "TF") %>%
  pivot_longer(
    cols = -TF,
    names_to = "comparison",
    values_to = "anova_FDR")

t2 <- t2 %>%
  rownames_to_column(var = "TF") %>%
  pivot_longer(
    cols = -TF,
    names_to = "comparison",
    values_to = "variance_explained")

t2$TF <- ifelse(t2$TF == "All 4\n(combinatorial)", "combined", t2$TF)

t2 <- t2 %>%
  left_join(t, by = c("TF", "comparison"))

write.csv(t2, file = "./analysis/csvs_for_supp_tables/hnDAM_TF_varpart_summary.csv", row.names = F)

##################################################

# BAM GEP

t <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/BAM_GEPs_df.rds")
write.csv(t, file = "./analysis/csvs_for_supp_tables/BAM_GEP.csv", row.names = F)

##################################################

# DAM replication datasets summary

t1 <- readRDS("./analysis/microglia/DAM_signature_discovery/PCA_GLM_summary.rds")
t2 <- readRDS("./analysis/microglia/DAM_signature_discovery/DAM_diffexp_GSEA_summary.rds")

t1$test <- gsub("group", "", t1$test)
t1$test <- rename_datasets(t1$test)

t2$test <- rename_datasets(t2$test)

summary <- t1 %>%
  left_join(t2, by = "test") %>%
  dplyr::select(-leadingEdge)


write.csv(summary, file = "./analysis/csvs_for_supp_tables/DAM_discovery_summary.csv", row.names = F)

##################################################

# mashr betas/lfsrs

# microglia

micro <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_mashr_object.rds")
micro_beta <- as.data.frame(micro$result$PosteriorMean)
micro_lfsr <- as.data.frame(micro$result$lfsr)

micro_lfsr <- micro_lfsr %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = -gene,
               names_to = "comparison",
               values_to = "lfsr")

micro_beta <- micro_beta %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = -gene,
               names_to = "comparison",
               values_to = "beta") %>%
  left_join(micro_lfsr, by = c("gene", "comparison"))

micro_beta$comparison <- rename_datasets(micro_beta$comparison)

write.csv(micro_beta, file = "./analysis/csvs_for_supp_tables/micro_mashr_results_all.csv", row.names = F)



# astrocytes

astro <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_diffexp_mashr_object.rds")
astro_beta <- as.data.frame(astro$result$PosteriorMean)
astro_lfsr <- as.data.frame(astro$result$lfsr)

astro_lfsr <- astro_lfsr %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = -gene,
               names_to = "comparison",
               values_to = "lfsr")

astro_beta <- astro_beta %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = -gene,
               names_to = "comparison",
               values_to = "beta") %>%
  left_join(astro_lfsr, by = c("gene", "comparison"))

astro_beta$comparison <- rename_datasets(astro_beta$comparison)

write.csv(astro_beta, file = "./analysis/csvs_for_supp_tables/astro_mashr_results_all.csv", row.names = F)



# oligodendrocytes

oligo <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_mashr_object.rds")
oligo_beta <- as.data.frame(oligo$result$PosteriorMean)
oligo_lfsr <- as.data.frame(oligo$result$lfsr)

oligo_lfsr <- oligo_lfsr %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = -gene,
               names_to = "comparison",
               values_to = "lfsr")

oligo_beta <- oligo_beta %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = -gene,
               names_to = "comparison",
               values_to = "beta") %>%
  left_join(oligo_lfsr, by = c("gene", "comparison"))

oligo_beta$comparison <- rename_datasets(oligo_beta$comparison)

write.csv(oligo_beta, file = "./analysis/csvs_for_supp_tables/oligo_mashr_results_all.csv", row.names = F)

##################################################

# k-means enrichments

# microglia 

t <- readRDS("./analysis/microglia/differential_expression/microglia_gobp_c1_c2.rds")
t2 <- readRDS("./analysis/microglia/differential_expression/microglia_gobp_c4_c6.rds")


t <- t[t$pvalue < 0.05,]
t2 <- t2[t2$pvalue < 0.05,]

t$kmeans_cluster <- "c1 + c2"
t2$kmeans_cluster <- "c4 + c6"

all <- rbind(t, t2)

write.csv(all, file = "./analysis/csvs_for_supp_tables/micro_kmeans_gobp.csv", row.names = F)



# astrocytes 

t <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_gobp_c2.rds")
t2 <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_gobp_c5.rds")


t <- t[t$pvalue < 0.05,]
t2 <- t2[t2$pvalue < 0.05,]

t$kmeans_cluster <- "c2"
t2$kmeans_cluster <- "c5"

all <- rbind(t, t2)

write.csv(all, file = "./analysis/csvs_for_supp_tables/astro_kmeans_gobp.csv", row.names = F)



# oligodendrocytes 

t <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_gobp_c2.rds")
t2 <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_gobp_c3.rds")


t <- t[t$pvalue < 0.05,]
t2 <- t2[t2$pvalue < 0.05,]

t$kmeans_cluster <- "c2"
t2$kmeans_cluster <- "c3"

all <- rbind(t, t2)

write.csv(all, file = "./analysis/csvs_for_supp_tables/oligo_kmeans_gobp.csv", row.names = F)

##################################################

# DESeq results from replication studies

dds_list <- list.files("./analysis/microglia/DAM_signature_replication/DESeq_results/", recursive = F)

list <- list()

for (file in dds_list){
  t <- readRDS(paste0("./analysis/microglia/DAM_signature_replication/DESeq_results/", file))
  name <- gsub("_deseq.rds", "", file)
  list[[name]] <- t
}

results_list <- list()

for (i in seq_along(list)) {
  dds <- list[[i]]
  dataset_label <- paste0("Dataset_", i)  # You can replace with a name vector if you have names
  
  cat("\n---", dataset_label, "---\n")
  rn <- resultsNames(dds)
  print(rn)
  
  results_list[[dataset_label]] <- list()  # Initialize sublist for this dataset
  
  repeat {
    res_name <- readline(prompt = paste("Enter result name to extract from", dataset_label, "(or 'q' to skip to next dataset): "))
    
    if (tolower(res_name) == "q") {
      cat("Moving to next dataset.\n")
      break
    }
    
    if (!(res_name %in% rn)) {
      cat("Invalid result name. Try again.\n")
      next
    }
    
    # Extract and save
    res_df <- as.data.frame(results(dds, name = res_name))
    results_list[[dataset_label]][[res_name]] <- res_df
    cat(paste("Added result:", res_name, "\n"))
    
    # Ask if user wants to extract more results from this dataset
    more <- readline(prompt = "Do you want to extract another result from this dataset? (y/n): ")
    if (tolower(more) != "y") {
      break
    }
  }
}

flat_results <- list()

for (dataset_name in names(results_list)) {
  for (res_name in names(results_list[[dataset_name]])) {
    combined_name <- paste(dataset_name, res_name, sep = "_")
    df <- results_list[[dataset_name]][[res_name]]
    flat_results[[combined_name]] <- na.omit(df)
  }
}

names(flat_results) <- c("Green_2024_PFC_Braak", "Green_2024_PFC_CERAD", "Green_2024_PFC_MCI", "Green_2024_PFC_AD_dementia",
                         "Kamath_2022_SN_PD/DLB", 
                         "Macnair_2025_CTX_GM_NAGM", "Macnair_2025_CTX_GM_GML",
                         "Macnair_2025_CTX_WM_NAWM", "Macnair_2025_CTX_WM_AL", "Macnair_2025_CTX_WM_CAL",
                         "Macnair_2025_CTX_WM_CIL", "Macnair_2025_CTX_WM_RL", 
                         "Martirosyan_2024_SN_PD",
                         "Mathys_2023_PFC_Braak", "Mathys_2023_PFC_CERAD", "Mathys_2023_PFC_MCI", "Mathys_2023_PFC_AD_dementia",
                         "Petrescu_2025_BA44_ALS", "Petrescu_2025_BA44_ALS_borderline", "Petrescu_2025_BA44_ALSci",
                         "Petrescu_2025_BA46_ALS", "Petrescu_2025_BA46_ALS_borderline", "Petrescu_2025_BA46_ALSci",
                         "Prudencio_2015_CB_sALS", "Prudencio_2015_CB_C9ALS", "Prudencio_2015_FC_sALS", "Prudencio_2015_FC_C9ALS",
                         "Rexach_2024_INS_AD", "Rexach_2024_INS_PSP", "Rexach_2024_INS_PiD",
                         "Gabitto_2024_PFC_Braak", "Gabitto_2024_PFC_CERAD", "Gabitto_2024_PFC_AD_dementia",
                         "Wang_2024_SN_PD/PDD",
                         "Zelic_2025_SpC_ALS")

for (test in names(flat_results)){
  df <- flat_results[[test]]
  test <- gsub("/", "_", test)
  write.csv(df, paste0("./analysis/csvs_for_supp_tables/DAM_replication_DESeq_results/", test, ".csv"), row.names = T)
}

##################################################

# num cells per subcluster per condition

donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")


# microglia

micro <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")

micro_cts <- micro %>%
  left_join(donor_meta, by = "donor")

micro_cts <- micro_cts %>%
  group_by(dataset, region, group, cluster_anno) %>%
  summarise(n = n())


# astrocytes

astro <- readRDS("./analysis/astrocytes/all_astros_metadata_ANNOTATED.rds")

astro_cts <- astro %>%
  left_join(donor_meta, by = "donor")

astro_cts <- astro_cts %>%
  group_by(dataset, region, group, cluster_anno) %>%
  summarise(n = n())


# oligodendrocytes

oligo <- readRDS("./analysis/oligodendrocytes/oligos_celllevel_metadata_ANNOTATED.rds")

oligo_cts <- oligo %>%
  left_join(donor_meta, by = "donor")

oligo_cts <- oligo_cts %>%
  group_by(dataset, region, group, cluster_anno) %>%
  summarise(n = n())


cts_merged <- rbind(micro_cts, astro_cts, oligo_cts)
cts_merged$dataset <- ifelse(cts_merged$dataset == "AMP-PD", "NM_2024",
                             ifelse(cts_merged$dataset == "Pineda", "Pineda_2024",
                                    ifelse(cts_merged$dataset == "Mathys", "Mathys_2024", "Gerrits_2022")))

cts_merged$region <- ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "FC/PFC", "PFC", 
                            ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "OC/V1", "PFC", 
                                   ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "TC/MTG", "MTG",
                                          ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "FC/PFC", "FC", 
                                                 ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "OC/V1", "OC", 
                                                        ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "TC/MTG", "TC", cts_merged$region))))))

colnames(cts_merged) <- c("dataset", "region", "disease_group", "cluster", "n_cells")

write.csv(cts_merged, file = "./analysis/csvs_for_supp_tables/cells_by_dataset_region_group.csv", row.names = F)

##################################################

# num cells per subcluster per condition, cortical/subcortical astros only

donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")


# cortical

cortical <- readRDS("./analysis/astrocytes/cortical_astros_metadata_ANNOTATED.rds")

cortical <- cortical %>%
  left_join(donor_meta, by = "donor")

cortical <- cortical %>%
  group_by(dataset, region, group, cluster_anno) %>%
  summarise(n = n())


# subcortical

subcortical <- readRDS("./analysis/astrocytes/subcortical_astros_metadata_ANNOTATED.rds")

subcortical <- subcortical %>%
  left_join(donor_meta, by = "donor")

subcortical <- subcortical %>%
  group_by(dataset, region, group, cluster_anno) %>%
  summarise(n = n())



cts_merged <- rbind(cortical, subcortical)
cts_merged$dataset <- ifelse(cts_merged$dataset == "AMP-PD", "NM_2024",
                             ifelse(cts_merged$dataset == "Pineda", "Pineda_2024",
                                    ifelse(cts_merged$dataset == "Mathys", "Mathys_2024", "Gerrits_2022")))

cts_merged$region <- ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "FC/PFC", "PFC", 
                            ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "OC/V1", "PFC", 
                                   ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "TC/MTG", "MTG",
                                          ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "FC/PFC", "FC", 
                                                 ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "OC/V1", "OC", 
                                                        ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "TC/MTG", "TC", cts_merged$region))))))

colnames(cts_merged) <- c("dataset", "region", "disease_group", "cluster", "n_cells")

write.csv(cts_merged, file = "./analysis/csvs_for_supp_tables/cells_by_dataset_region_group_astros_subclustered.csv", row.names = F)

##################################################

# overall discovery dataset summary

# microglia

donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")


# microglia

micro <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")

micro <- micro %>%
  left_join(donor_meta, by = "donor")

micro_cts <- micro %>%
  group_by(dataset, region, group) %>%
  summarise(n_cells = n(),
            n_donors = n_distinct(donor))

micro_cts$celltype <- "Microglia"


# astrocytes

astro <- readRDS("./analysis/astrocytes/all_astros_metadata_ANNOTATED.rds")

astro <- astro %>%
  left_join(donor_meta, by = "donor")

astro_cts <- astro %>%
  group_by(dataset, region, group) %>%
  summarise(n_cells = n(),
            n_donors = n_distinct(donor))

astro_cts$celltype <- "Astrocytes"


# oligodendrocytes

oligo <- readRDS("./analysis/oligodendrocytes/oligos_celllevel_metadata_ANNOTATED.rds")

oligo <- oligo %>%
  left_join(donor_meta, by = "donor")

oligo_cts <- oligo %>%
  group_by(dataset, region, group) %>%
  summarise(n_cells = n(),
            n_donors = n_distinct(donor))

oligo_cts$celltype <- "Oligodendrocytes"


cts_merged <- rbind(micro_cts, astro_cts, oligo_cts)
cts_merged$dataset <- ifelse(cts_merged$dataset == "AMP-PD", "NM_2024",
                             ifelse(cts_merged$dataset == "Pineda", "Pineda_2024",
                                    ifelse(cts_merged$dataset == "Mathys", "Mathys_2024", "Gerrits_2022")))

cts_merged$region <- ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "FC/PFC", "PFC", 
                            ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "OC/V1", "PFC", 
                                   ifelse((!cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "TC/MTG", "MTG",
                                          ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "FC/PFC", "FC", 
                                                 ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "OC/V1", "OC", 
                                                        ifelse((cts_merged$dataset == "Gerrits_2022") & cts_merged$region == "TC/MTG", "TC", cts_merged$region))))))

colnames(cts_merged) <- c("dataset", "region", "disease_group", "n_cells", "n_unique_donors", "celltype")

write.csv(cts_merged, file = "./analysis/csvs_for_supp_tables/cells_donors_by_dataset_region_group.csv", row.names = F)

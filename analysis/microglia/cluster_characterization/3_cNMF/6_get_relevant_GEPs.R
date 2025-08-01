library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################

Gerrits_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_gene_scores.rds")
Mathys_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_gene_scores.rds")
Pineda_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_gene_scores.rds")

##################################################
##################################################
##################################################

# getting consensus DAM signature

##################################################

# genes expressed in all 3 datasets
genes <- intersect(rownames(Gerrits_allgenes), intersect(rownames(Mathys_allgenes), rownames(Pineda_allgenes)))


# get DAM signatures w/ enrichments from each dataset
Gerrits_GEP5 <- Gerrits_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP5) %>%
  dplyr::rename(Gerrits_GEP5 = GEP5)

Mathys_GEP3 <- Mathys_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP3) %>%
  dplyr::rename(Mathys_GEP3 = GEP3)

Pineda_GEP5 <- Pineda_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP5) %>%
  dplyr::rename(Pineda_GEP5 = GEP5)

Gerrits_GEP5 <- Gerrits_GEP5[Gerrits_GEP5$gene %in% genes,]
Mathys_GEP3 <- Mathys_GEP3[Mathys_GEP3$gene %in% genes,]
Pineda_GEP5 <- Pineda_GEP5[Pineda_GEP5$gene %in% genes,]


# join
DAM_GEPs <- Gerrits_GEP5 %>%
  left_join(Mathys_GEP3, by = "gene") %>%
  left_join(Pineda_GEP5, by = "gene")

##################################################
##################################################
##################################################

# genes w/ average z-scored enrichment across studies > 5 

DAM_GEPs$Gerrits_z <- scale(DAM_GEPs$Gerrits_GEP5)
DAM_GEPs$Mathys_z <- scale(DAM_GEPs$Mathys_GEP3)
DAM_GEPs$Pineda_z <- scale(DAM_GEPs$Pineda_GEP5)

DAM_GEPs$avg_z <- rowMeans(DAM_GEPs[, c("Gerrits_z", "Mathys_z", "Pineda_z")])

DAM_GEPs <- DAM_GEPs %>%
  arrange(desc(avg_z))

DAM_genes_zscore <- DAM_GEPs$gene[DAM_GEPs$avg_z >= 5]

saveRDS(DAM_genes_zscore, file = "./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
saveRDS(DAM_GEPs, file = "./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_GEPs_zscore_df.rds")
writeLines(DAM_genes_zscore, "./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.txt")


# write outputs for homer motif enrichment analysis
background <- DAM_GEPs$gene

writeLines(DAM_genes_zscore, "./analysis/microglia/cluster_characterization/homer_DAM_sig/DAM_genes.txt")
writeLines(background, "./analysis/microglia/cluster_characterization/homer_DAM_sig/background_genes.txt")

##################################################
##################################################
##################################################

# getting consensus homeostatic microglia signature

##################################################

Gerrits_GEP1 <- Gerrits_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP1) %>%
  dplyr::rename(Gerrits_GEP1 = GEP1)

Mathys_GEP1 <- Mathys_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP1) %>%
  dplyr::rename(Mathys_GEP1 = GEP1)

Pineda_GEP1 <- Pineda_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP1) %>%
  dplyr::rename(Pineda_GEP1 = GEP1)

Gerrits_GEP1 <- Gerrits_GEP1[Gerrits_GEP1$gene %in% genes,]
Mathys_GEP1 <- Mathys_GEP1[Mathys_GEP1$gene %in% genes,]
Pineda_GEP1 <- Pineda_GEP1[Pineda_GEP1$gene %in% genes,]


# join
homeo_GEPs <- Gerrits_GEP1 %>%
  left_join(Mathys_GEP1, by = "gene") %>%
  left_join(Pineda_GEP1, by = "gene")

##################################################

# genes w/ average z-scored enrichment across studies > 5 

homeo_GEPs$Gerrits_z <- scale(homeo_GEPs$Gerrits_GEP1)
homeo_GEPs$Mathys_z <- scale(homeo_GEPs$Mathys_GEP1)
homeo_GEPs$Pineda_z <- scale(homeo_GEPs$Pineda_GEP1)

homeo_GEPs$avg_z <- rowMeans(homeo_GEPs[, c("Gerrits_z", "Mathys_z", "Pineda_z")])

homeo_GEPs <- homeo_GEPs %>%
  arrange(desc(avg_z))

homeo_genes_zscore <- homeo_GEPs$gene[homeo_GEPs$avg_z >= 5]

saveRDS(homeo_genes_zscore, file = "./analysis/microglia/cluster_characterization/final_gene_signatures/homeo_genes_zscore.rds")
saveRDS(homeo_GEPs, file = "./analysis/microglia/cluster_characterization/final_gene_signatures/homeo_GEPs_zscore_df.rds")
writeLines(homeo_genes_zscore, "./analysis/microglia/cluster_characterization/final_gene_signatures/homeo_genes_zscore.txt")

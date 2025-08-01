library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################
##################################################
##################################################

Gerrits_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_gene_scores.rds")
Mathys_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_gene_scores.rds")
Pineda_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_gene_scores.rds")
AMPPD_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_gene_scores.rds")

#################################################
##################################################
##################################################

# getting consensus BAM signature

##################################################

# genes expressed in all 3 datasets
genes <- intersect(rownames(Gerrits_allgenes), intersect(rownames(Mathys_allgenes), intersect(rownames(Pineda_allgenes), rownames(AMPPD_allgenes))))


# get BAM signatures w/ enrichments from each dataset
Gerrits_GEP9 <- Gerrits_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP9) %>%
  dplyr::rename(Gerrits_GEP9 = GEP9)

Mathys_GEP6 <- Mathys_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP6) %>%
  dplyr::rename(Mathys_GEP6 = GEP6)

Pineda_GEP6 <- Pineda_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP6) %>%
  dplyr::rename(Pineda_GEP6 = GEP6)

AMPPD_GEP4 <- AMPPD_allgenes %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, GEP4) %>%
  dplyr::rename(AMPPD_GEP4 = GEP4)

Gerrits_GEP9 <- Gerrits_GEP9[Gerrits_GEP9$gene %in% genes,]
Mathys_GEP6 <- Mathys_GEP6[Mathys_GEP6$gene %in% genes,]
Pineda_GEP6 <- Pineda_GEP6[Pineda_GEP6$gene %in% genes,]
AMPPD_GEP4 <- AMPPD_GEP4[AMPPD_GEP4$gene %in% genes,]

# join
BAM_GEPs <- Gerrits_GEP9 %>%
  left_join(Mathys_GEP6, by = "gene") %>%
  left_join(Pineda_GEP6, by = "gene") %>%
  left_join(AMPPD_GEP4, by = "gene")

##################################################
##################################################
##################################################

# genes w/ average z-scored enrichment across studies > 5 

BAM_GEPs$Gerrits_z <- scale(BAM_GEPs$Gerrits_GEP9)
BAM_GEPs$Mathys_z <- scale(BAM_GEPs$Mathys_GEP6)
BAM_GEPs$Pineda_z <- scale(BAM_GEPs$Pineda_GEP6)
BAM_GEPs$AMPPD_z <- scale(BAM_GEPs$AMPPD_GEP4)

BAM_GEPs$avg_z <- rowMeans(BAM_GEPs[, c("Gerrits_z", "Mathys_z", "Pineda_z", "AMPPD_z")])

BAM_GEPs <- BAM_GEPs %>%
  arrange(desc(avg_z))

BAM_genes_zscore <- BAM_GEPs$gene[BAM_GEPs$avg_z >= 5]

saveRDS(BAM_genes_zscore, file = "./analysis/microglia/cluster_characterization/final_gene_signatures/BAM_genes_zscore.rds")
saveRDS(BAM_GEPs, file = "./analysis/microglia/cluster_characterization/final_gene_signatures/BAM_GEPs_df.rds")

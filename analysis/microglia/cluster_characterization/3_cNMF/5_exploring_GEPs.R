library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(Seurat)
library(scCustomize)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################

micro <- readRDS("./analysis/microglia/seurat_objects/microglia_annotated.rds")

AMPPD_GEPs <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_top200_list.rds")
Gerrits_GEPs <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_top200_list.rds")
Mathys_GEPs <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_top200_list.rds")
Pineda_GEPs <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_top200_list.rds")

AMPPD_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_gene_scores.rds")
Gerrits_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_gene_scores.rds")
Mathys_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_gene_scores.rds")
Pineda_allgenes <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_gene_scores.rds")

AMPPD_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_avg_cluster_usage.rds")
Gerrits_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_avg_cluster_usage.rds")
Mathys_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_avg_cluster_usage.rds")
Pineda_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_avg_cluster_usage.rds")

##################################################

col_fun = colorRamp2(breaks = c(0, 1), colors = c("white", "red"))

Heatmap(as.matrix(AMPPD_usages),
        cluster_rows = F,
        cluster_columns = F,
        col = col_fun,
        name = "GEP usage")

# Gerrits -- GEP1 = homeostatic, GEP5 = DAM
# Mathys -- GEP1 = homeostatic, GEP3 = DAM
# Pineda -- GEP1 = homeostatic, GEP5 = DAM

##################################################

Gerrits_DAM_genes <- Gerrits_GEPs$GEP5$gene#[1:100]
Mathys_DAM_genes <- Mathys_GEPs$GEP3$gene#[1:100]
Pineda_DAM_genes <- Pineda_GEPs$GEP5$gene#[1:100]

DAM_genes <- intersect(Gerrits_DAM_genes, intersect(Mathys_DAM_genes, Pineda_DAM_genes))
#writeLines(DAM_genes, "./analysis/microglia/misc_files/DAM_genes.txt")


gobp <- enrichGO(gene = DAM_genes, 
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP", 
                 keyType = "SYMBOL", 
                 readable = TRUE,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05)
gobp_df <- gobp_df@result
gobp_df <- gobp_df[gobp_df$qvalue < 0.05, ]

mutate(gobp, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

dotplot(gobp, showCategory = 20)

##################################################

# determine genes for GEPs relatively rather than hard 1:200
Gerrits_DAM_genes_relative <- rownames(Gerrits_allgenes)[Gerrits_allgenes$GEP5 >= quantile(Gerrits_allgenes$GEP5, probs = 0.99)]
Mathys_DAM_genes_relative <- rownames(Mathys_allgenes)[Mathys_allgenes$GEP3 >= quantile(Mathys_allgenes$GEP3, probs = 0.99)]
Pineda_DAM_genes_relative <- rownames(Pineda_allgenes)[Pineda_allgenes$GEP5 >= quantile(Pineda_allgenes$GEP5, probs = 0.99)]

DAM_genes_relative <- intersect(Gerrits_DAM_genes_relative, intersect(Mathys_DAM_genes_relative, Pineda_DAM_genes_relative))


gobp_relative <- enrichGO(gene = DAM_genes, 
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          keyType = "SYMBOL",
                          readable = TRUE,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)
gobp_relative <- gobp_relative@result
gobp_relative <- gobp_relative[gobp_relative$qvalue < 0.05, ]

##################################################

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


DAM_GEPs_merged <- Mathys_GEP3 %>%
  left_join(Gerrits_GEP5, by = "gene") %>%
  left_join(Pineda_GEP5, by = "gene")


# rank genes w/in GEPs and plot concordance of rankings
DAM_GEPs_merged_filtered <- na.omit(DAM_GEPs_merged)

DAM_GEPs_ranked <- DAM_GEPs_merged_filtered
DAM_GEPs_ranked[, -1] <- apply(DAM_GEPs_ranked[, -1], 2, function(x) rank(-x, ties.method = "min"))
DAM_GEPs_ranked <- DAM_GEPs_ranked %>%
  dplyr::arrange(Mathys_GEP3)


ggplot(DAM_GEPs_ranked[(nrow(DAM_GEPs_ranked)-100):nrow(DAM_GEPs_ranked),], aes(x = Mathys_GEP3, y = Pineda_GEP5)) + 
  geom_point()
ggplot(DAM_GEPs_ranked[1:200,], aes(x = Mathys_GEP3, y = Pineda_GEP5)) + 
  geom_point()




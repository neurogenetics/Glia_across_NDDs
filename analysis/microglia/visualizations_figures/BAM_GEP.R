library(tidyverse)
library(Seurat)
library(scCustomize)
library(ggplot2)
library(clusterProfiler)
library(rrvgo)
library(org.Hs.eg.db)
library(svglite)
library(ComplexHeatmap)
library(circlize)
library(Cairo)
library(eulerr)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################

source("./code_organized/functions/plotting_colors.R")

##################################################
##################################################
##################################################

# plotting z-score enrichments of individual genes 

BAM_genes_df <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/BAM_GEPs_df.rds")
DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")

BAM_genes_df <- BAM_genes_df %>%
  arrange(avg_z)

p1 <- ggplot(BAM_genes_df, aes(x = 1:nrow(BAM_genes_df), y = avg_z, color = avg_z > 5)) +
  geom_point() +
  labs(y = "Average cross-study\nz-score gene enrichment",
       x = "Gene rank") +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "purple", linewidth = 1.5) + 
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "purple"))

ggsave(p1, filename = "./analysis/microglia/plots_final/BAM_genes_zscore_enrichment.png", width = 4.5, height = 3.2, dpi = 600)

##################################################

# overlay enrichment of hnDAM genes onto BAM ranks

BAM_genes_df$rank <- 1:nrow(BAM_genes_df)

p3 <- ggplot() +
  geom_point(data = subset(BAM_genes_df, !(gene %in% DAM_genes)), 
             aes(x = rank, y = avg_z), color = "gray") +
  geom_point(data = subset(BAM_genes_df, gene %in% DAM_genes), 
             aes(x = rank, y = avg_z), color = "yellow2") +
  labs(y = "Average cross-study\nz-score gene enrichment",
       x = "BAM gene rank") +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "purple", linewidth = 1.5)

ggsave(p3, filename = "./analysis/microglia/plots_final/BAM_genes_zscore_enrichment_hnDAM_overlaid.png", width = 4.5, height = 3.2, dpi = 600)

##################################################
##################################################
##################################################

# heatmap of gene expression of BAM genes across clusters

BAM_genes_up <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/BAM_genes_zscore.rds")

micro <- readRDS("./analysis/microglia/seurat_objects/microglia_annotated.rds")
micro[["RNA"]]$data <- as(object = micro[["RNA"]]$data, Class = "dgCMatrix")

avgexp <- as.data.frame(AverageExpression(micro, assays = "RNA", group.by = "cluster_anno", 
                                          features = BAM_genes_up, layer = "data"))

saveRDS(avgexp, file = "./analysis/microglia/DAM_signature_discovery/BAM_genes_avgexp_across_clusters.rds")

#########################

avgexp <- readRDS("./analysis/microglia/DAM_signature_discovery/BAM_genes_avgexp_across_clusters.rds")

expr_scaled <- as.data.frame(t(apply(avgexp, 1, function(x) scale(x))))
colnames(expr_scaled) <- colnames(avgexp)
colnames(expr_scaled) <- gsub("RNA.", "", colnames(expr_scaled))
colnames(expr_scaled) <- gsub("\\.", "_", colnames(expr_scaled))

expr_scaled <- expr_scaled[, c("Micro_Homeo", 
                               "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                               "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                               "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN",
                               "BAM", "Monocyte", "Lymphocyte")]


col_fun = colorRamp2(breaks = c(-2, 0, 4), colors = c("#05409e", "white", "#e64a02"))

genes_label = c(BAM_genes_up[1:10], "LYVE1")

col_idx <- which(rownames(expr_scaled) %in% genes_label)

col_anno <- columnAnnotation(
  mark = anno_mark(
    at = col_idx,
    side = "bottom",
    labels = rownames(expr_scaled)[col_idx]))

ht = Heatmap(as.matrix(t(expr_scaled)),
             cluster_rows = F,
             cluster_columns = T,
             col = col_fun,
             name = "Scaled expression",
             show_column_dend = F,
             bottom_annotation = col_anno,
             show_row_names = F,
             show_column_names = F,
             heatmap_legend_param = list(title_position = "topleft",
                                         legend_direction = "horizontal",
                                         legend_width = unit(3, "cm"),
                                         title_gp = gpar(fontsize = 12)))
draw(ht, heatmap_legend_side = "bottom")


svglite(filename = "./analysis/microglia/plots_final/BAM_sig_genes_expr.svg", width = 9.6, height = 4.3)
draw(ht)
dev.off()

##################################################
##################################################
##################################################

# show relevant GEP usages from cNMF

AMPPD_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_avg_cluster_usage.rds")
Gerrits_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_avg_cluster_usage.rds")
Mathys_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_avg_cluster_usage.rds")
Pineda_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_avg_cluster_usage.rds")

AMPPD_usages <- AMPPD_usages %>%
  rownames_to_column(var = "celltype") %>%
  dplyr::select(celltype, GEP4) %>%
  dplyr::rename(NM_2024_GEP4 = GEP4)

Gerrits_usages <- Gerrits_usages %>%
  rownames_to_column(var = "celltype") %>%
  dplyr::select(celltype, GEP9) %>%
  dplyr::rename(Gerrits_2022_GEP9 = GEP9)

Mathys_usages <- Mathys_usages %>%
  rownames_to_column(var = "celltype") %>%
  dplyr::select(celltype, GEP6) %>%
  dplyr::rename(Mathys_2024_GEP6 = GEP6)

Pineda_usages <- Pineda_usages %>%
  rownames_to_column(var = "celltype") %>%
  dplyr::select(celltype, GEP6) %>%
  dplyr::rename(Pineda_2024_GEP6 = GEP6)


usages_merged <- AMPPD_usages %>%
  left_join(Mathys_usages, by = "celltype") %>%
  left_join(Pineda_usages, by = "celltype") %>%
  left_join(Gerrits_usages, by = "celltype") %>%
  column_to_rownames(var = "celltype")


cluster_order <- c("Micro_Homeo", 
                   "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                   "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                   "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN",
                   "BAM", "Monocyte", "Lymphocyte")

usages_merged <- usages_merged[cluster_order, ]

usages_merged <- usages_merged[c("Gerrits_2022_GEP9", "Mathys_2024_GEP6", "Pineda_2024_GEP6", "NM_2024_GEP4")]

col_fun = colorRamp2(breaks = c(0, 1), colors = c("white", "red"))

ht = Heatmap(as.matrix(usages_merged),
             cluster_rows = F,
             cluster_columns = F,
             col = col_fun,
             name = "GEP usage", 
             rect_gp = gpar(col = "black"),
             column_names_rot = 45,
             column_gap = unit(5, "mm"),
             column_title = NULL,
             column_names_side = "top", 
             heatmap_legend_param = list(title_position = "topleft",
                                         legend_direction = "horizontal",
                                         legend_width = unit(3, "cm"),
                                         title_gp = gpar(fontsize = 12)))

cairo_pdf("./analysis/microglia/plots_final/BAM_cNMF_GEP_usages.pdf", width = 4.5, height = 5.1, family = "Arial")
draw(ht, heatmap_legend_side = "bottom")
dev.off()

##################################################
##################################################
##################################################

# seurat module score 

BAM_genes_up <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/BAM_genes_zscore.rds")

micro <- readRDS("./analysis/microglia/seurat_objects/microglia_annotated.rds")

cluster_order <- c("Micro_Homeo", 
                   "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                   "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                   "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN",
                   "BAM", "Monocyte", "Lymphocyte")

module <- list(BAM_genes_up)

micro <- AddModuleScore(micro, features = module)


micro$cluster_anno <- factor(micro$cluster_anno, levels = cluster_order)
Idents(micro) <- "cluster_anno"

p2 <- VlnPlot_scCustom(micro, group.by = "cluster_anno", features = "Cluster1", 
                      pt.size = 0, colors_use = microglia_colors) + 
  labs(y = "BAM signature\nmodule enrichment",
       title = NULL) +
  geom_boxplot(outliers = F, width = 0.15, fill = "lightgrey") + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12)) +
  NoLegend()

svglite("./analysis/microglia/plots_final/BAM_seurat_modulescore_vlnplot.svg", width = 8.3, height = 4)
plot(p2)
dev.off()

##################################################
##################################################
##################################################

# overlap b/w DAM and BAM

BAM_genes_up <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/BAM_genes_zscore.rds")
DAM_genes_up <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")

intersect(BAM_genes_up, DAM_genes_up)

plot(venn(list(BAM = BAM_genes_up, hnDAM = DAM_genes_up)))

png("./analysis/microglia/plots_final/DAM_BAM_overlap_venn.png", width = 4.8, height = 2.75, units = "in", res = 600)
plot(venn(list(BAM = BAM_genes_up, hnDAM = DAM_genes_up)), fills = list(fill = c("#968863", "yellow"), alpha = 0.6),
     labels = T, edges = T, quantities = T)
dev.off()


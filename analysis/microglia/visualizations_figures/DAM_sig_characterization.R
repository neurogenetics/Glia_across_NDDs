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

setwd("/data/ADRD/glia_across_NDDs/")

##################################################

source("./code_organized/functions/plotting_colors.R")

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
  dplyr::select(celltype, GEP1) %>%
  dplyr::rename(NM_2024_GEP1 = GEP1)

Gerrits_usages <- Gerrits_usages %>%
  rownames_to_column(var = "celltype") %>%
  dplyr::select(celltype, GEP5) %>%
  dplyr::rename(Gerrits_2022_GEP5 = GEP5)

Mathys_usages <- Mathys_usages %>%
  rownames_to_column(var = "celltype") %>%
  dplyr::select(celltype, GEP3) %>%
  dplyr::rename(Mathys_2024_GEP3 = GEP3)

Pineda_usages <- Pineda_usages %>%
  rownames_to_column(var = "celltype") %>%
  dplyr::select(celltype, GEP5) %>%
  dplyr::rename(Pineda_2024_GEP5 = GEP5)


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

usages_merged <- usages_merged[c("Gerrits_2022_GEP5", "Mathys_2024_GEP3", "Pineda_2024_GEP5", "NM_2024_GEP1")]

column_split <- ifelse(colnames(usages_merged) == "NM_2024_GEP1", "Group2", "Group1")

col_fun = colorRamp2(breaks = c(0, 1), colors = c("white", "red"))

ht = Heatmap(as.matrix(usages_merged),
        cluster_rows = F,
        cluster_columns = F,
        col = col_fun,
        name = "GEP usage", 
        rect_gp = gpar(col = "black"),
        column_names_rot = 45,
        column_split = column_split,
        column_gap = unit(5, "mm"),
        column_title = NULL,
        column_names_side = "top", 
        heatmap_legend_param = list(title_position = "topleft",
                                    legend_direction = "horizontal",
                                    legend_width = unit(3, "cm"),
                                    title_gp = gpar(fontsize = 12)))

cairo_pdf("./analysis/microglia/plots_final/DAM_cNMF_GEP_usages.pdf", width = 4.8, height = 5.1, family = "Arial")
draw(ht, heatmap_legend_side = "bottom")
dev.off()

##################################################
##################################################
##################################################

# seurat module score 

DAM_genes_up <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")

micro <- readRDS("./analysis/microglia/seurat_objects/microglia_annotated.rds")

cluster_order <- c("Micro_Homeo", 
                   "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                   "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                   "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN",
                   "BAM", "Monocyte", "Lymphocyte")

module <- list(DAM_genes_up)

micro <- AddModuleScore(micro, features = module)


micro$cluster_anno <- factor(micro$cluster_anno, levels = cluster_order)
Idents(micro) <- "cluster_anno"

p <- VlnPlot_scCustom(micro, group.by = "cluster_anno", features = "Cluster1", 
                      pt.size = 0, colors_use = microglia_colors) + 
  labs(y = "DAM signature\nmodule enrichment",
       title = NULL) +
  geom_boxplot(outliers = F, width = 0.15, fill = "lightgrey") + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12)) +
  NoLegend()

svglite("./analysis/microglia/plots_final/DAM_seurat_modulescore_vlnplot.svg", width = 8.3, height = 4)
plot(p)
dev.off()

##################################################
##################################################
##################################################

# gene ontologies

gobp_DAM <- enrichGO(gene = DAM_genes_up,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     keyType = "SYMBOL",
                     readable = TRUE,
                     pAdjustMethod = "BH")
gobp_DAM <- gobp_DAM@result
saveRDS(gobp_DAM, file = "./analysis/microglia/cluster_characterization/cNMF_by_dataset/DAM_sig_GOBP_all.rds")


gobp_DAM <- gobp_DAM[gobp_DAM$p.adjust < 0.05, ]


simMatrix <- calculateSimMatrix(gobp_DAM$ID, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
scores <- setNames(-log10(gobp_DAM$qvalue), gobp_DAM$ID)
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")

png(filename = "./analysis/microglia/plots_final/DAM_module_GO.png", width = 4, height = 3, units = "in", res = 600)
treemapPlot(reducedTerms)
dev.off()

##################################################
##################################################
##################################################

# plotting z-score enrichments of individual genes 

DAM_genes_df <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_GEPs_zscore_df.rds")

DAM_genes_df <- DAM_genes_df %>%
  arrange(avg_z)

p <- ggplot(DAM_genes_df, aes(x = 1:nrow(DAM_genes_df), y = avg_z, color = avg_z > 5)) +
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

svglite("./analysis/microglia/plots_final/DAM_genes_zscore_enrichment.svg", width = 4, height = 2.25)
plot(p)
dev.off()

##################################################
##################################################
##################################################

# heatmap of gene expression of DAM genes across clusters

DAM_genes_up <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")

micro <- readRDS("./analysis/microglia/seurat_objects/microglia_annotated.rds")
micro[["RNA"]]$data <- as(object = micro[["RNA"]]$data, Class = "dgCMatrix")

avgexp <- as.data.frame(AverageExpression(micro, assays = "RNA", group.by = "cluster_anno", 
                                          features = DAM_genes_up, layer = "data"))

saveRDS(avgexp, file = "./analysis/microglia/DAM_signature_discovery/DAM_genes_avgexp_across_clusters.rds")

#########################

avgexp <- readRDS("./analysis/microglia/DAM_signature_discovery/DAM_genes_avgexp_across_clusters.rds")

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

genes_label = c(DAM_genes_up[1:10], "SPP1")

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
draw(ht)


svglite(filename = "./analysis/microglia/plots_final/DAM_sig_genes_expr.svg", width = 9.6, height = 4.3)
draw(ht)
dev.off()

##################################################
##################################################
##################################################

# show all GEP usages from cNMF

AMPPD_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/AMPPD_microglia_GEPs_avg_cluster_usage.rds")
Gerrits_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Gerrits_microglia_GEPs_avg_cluster_usage.rds")
Mathys_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Mathys_microglia_GEPs_avg_cluster_usage.rds")
Pineda_usages <- readRDS("./analysis/microglia/cluster_characterization/cNMF_by_dataset/final_cNMF_GEPs/Pineda_microglia_GEPs_avg_cluster_usage.rds")


cluster_order <- c("Micro_Homeo", 
                   "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                   "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                   "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN",
                   "BAM", "Monocyte", "Lymphocyte")

AMPPD_usages <- AMPPD_usages[cluster_order, ]
Gerrits_usages <- Gerrits_usages[cluster_order, ]
Mathys_usages <- Mathys_usages[cluster_order, ]
Pineda_usages <- Pineda_usages[cluster_order, ]

Gerrits_usages <- na.omit(Gerrits_usages)
Mathys_usages <- na.omit(Mathys_usages)
Pineda_usages <- na.omit(Pineda_usages)


col_fun = colorRamp2(breaks = c(0, 1), colors = c("white", "red"))

ht1 = Heatmap(as.matrix(AMPPD_usages),
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

ht2 = Heatmap(as.matrix(Gerrits_usages),
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

ht3 = Heatmap(as.matrix(Mathys_usages),
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

ht4 = Heatmap(as.matrix(Pineda_usages),
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


cairo_pdf("./analysis/microglia/plots_final/AMPPD_all_GEP_usages.pdf", width = 10.4, height = 5.1, family = "Arial")
draw(ht1, heatmap_legend_side = "bottom")
dev.off()

cairo_pdf("./analysis/microglia/plots_final/Gerrits_all_GEP_usages.pdf", width = 10.4, height = 5.1, family = "Arial")
draw(ht2, heatmap_legend_side = "bottom")
dev.off()

cairo_pdf("./analysis/microglia/plots_final/Mathys_all_GEP_usages.pdf", width = 10.4, height = 5.1, family = "Arial")
draw(ht3, heatmap_legend_side = "bottom")
dev.off()

cairo_pdf("./analysis/microglia/plots_final/Pineda_all_GEP_usages.pdf", width = 10.4, height = 5.1, family = "Arial")
draw(ht4, heatmap_legend_side = "bottom")
dev.off()

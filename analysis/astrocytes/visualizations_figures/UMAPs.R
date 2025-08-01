library(Seurat)
library(scCustomize)
library(tidyverse)
library(cowplot)

setwd("/data/ADRD/glia_across_NDDs/")

########################################
########################################
########################################

cortical <- readRDS("./analysis/astrocytes/seurat_objects/cortical_astrocytes_clustered_ANNOTATED.rds")
subcortical <- readRDS("./analysis/astrocytes/seurat_objects/subcortical_astrocytes_clustered_ANNOTATED.rds")
astros_all <- readRDS("./analysis/astrocytes/seurat_objects/all_astrocytes_clustered_ANNOTATED.rds")

source("./code_organized/functions/plotting_colors.R")

########################################
########################################
########################################

# overall UMAP of all astros

Idents(astros_all) <- "cluster_anno"

cluster_order <- rev(c("Astro_Protoplasmic_GRM3", "Astro_Fibrous_DLCK1", "Astro_Fibrous_GRIA1",
                       "Astro_KCND2", "Astro_Reactive_SERPINA3", "Astro_LAMA2", "Astro_Stress_HSPH1"))

Idents(astros_all) <- factor(Idents(astros_all), levels = rev(cluster_order))

p1 <- DimPlot_scCustom(astros_all, repel = T, label.box = T, raster = F, colors_use = astrocyte_colors) +
  theme(text = element_text(family = "Arial"), 
        legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  coord_fixed()

ggsave(plot = p1, filename = "./analysis/astrocytes/plots_final/all_astros_umap_clusters.png", width = 9, height = 6.5, dpi = 600)

########################################
########################################
########################################

# all astros, highlighting cortical/subcortical

astros_all$broad_region <- ifelse(astros_all$region %in% c("DMV", "GPi", "TH"), "Subcortical", "Cortical")

astros_all$broad_region <- factor(astros_all$broad_region, levels = c("Cortical", "Subcortical"))

DimPlot_scCustom(astros_all, group.by = "broad_region", colors_use = c("#9e1f63", "#6fc251"), raster = F) + coord_fixed()

Idents(astros_all) <- "broad_region"

p2 <- Cluster_Highlight_Plot(astros_all, cluster_name = "Cortical", highlight_color = "#9e1f63", raster = F) + coord_fixed() + NoLegend() +
  ggtitle("Cortical") + 
  theme(text = element_text(family = "Arial"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 25))
p3 <- Cluster_Highlight_Plot(astros_all, cluster_name = "Subcortical", highlight_color = "#6fc251", raster = F) + coord_fixed() + NoLegend() +
  ggtitle("Subcortical") + 
  theme(text = element_text(family = "Arial"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 25))


ggsave(p2, filename = "./analysis/astrocytes/plots_final/astros_highlight_cortical.png", height = 7, width = 8.5, dpi = 600)
ggsave(p3, filename = "./analysis/astrocytes/plots_final/astros_highlight_subcortical.png", height = 7, width = 8.5, dpi = 600)

########################################
########################################
########################################

# cortical astros clusters

p4 <- DimPlot_scCustom(cortical, repel = F, label.box = T, raster = F, colors_use = astrocyte_colors) +
  theme(text = element_text(family = "Arial"), 
        legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  coord_fixed()

ggsave(plot = p4, filename = "./analysis/astrocytes/plots_final/cortical_astros_umap_clusters.png", width = 9, height = 7.5, dpi = 600)

########################################
########################################
########################################

# subcortical astros clusters

p5 <- DimPlot_scCustom(subcortical, repel = T, label.box = T, raster = F, colors_use = astrocyte_colors) +
  theme(text = element_text(family = "Arial"), 
        legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  coord_fixed()

ggsave(plot = p5, filename = "./analysis/astrocytes/plots_final/subcortical_astros_umap_clusters.png", width = 9, height = 7.5, dpi = 600)

########################################
########################################
########################################

# cortical astros, highlighting neocortex/allocortex

cortical$broad_region <- ifelse(cortical$region %in% c("EC", "HIP"), "Allocortical", "Neocortical")

cortical$broad_region <- factor(cortical$broad_region, levels = c("Neocortical", "Allocortical"))

DimPlot_scCustom(cortical, group.by = "broad_region", colors_use = c("#9e1f63", "#6fc251"), raster = F) + coord_fixed()

Idents(cortical) <- "broad_region"

p6 <- Cluster_Highlight_Plot(cortical, cluster_name = "Neocortical", highlight_color = "#4992eb", raster = F) + coord_fixed() + NoLegend() +
  ggtitle("Neocortical") + 
  theme(text = element_text(family = "Arial"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 25))
p7 <- Cluster_Highlight_Plot(cortical, cluster_name = "Allocortical", highlight_color = "#eb5449", raster = F) + coord_fixed() + NoLegend() +
  ggtitle("Allocortical") + 
  theme(text = element_text(family = "Arial"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 25))

p8 <- p6 / p7

ggsave(plot = p8, filename = "./analysis/astrocytes/plots_final/cortical_astros_neo_allo_highlight.png", width = 10, height = 15, dpi = 600)

########################################
########################################
########################################

# subcortical astros, highlighting regions

Idents(subcortical) <- "region"

p9 <- Cluster_Highlight_Plot(subcortical, cluster_name = "DMV", highlight_color = brain_region_colors, raster = F) + coord_fixed() + NoLegend() +
  ggtitle("DMV") + 
  theme(text = element_text(family = "Arial"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 25))
p10 <- Cluster_Highlight_Plot(subcortical, cluster_name = "GPi", highlight_color = brain_region_colors, raster = F) + coord_fixed() + NoLegend() +
  ggtitle("GPi") + 
  theme(text = element_text(family = "Arial"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 25))
p11 <- Cluster_Highlight_Plot(subcortical, cluster_name = "TH", highlight_color = brain_region_colors, raster = F) + coord_fixed() + NoLegend() +
  ggtitle("TH") + 
  theme(text = element_text(family = "Arial"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 25))

p12 <- p9/p10/p11

ggsave(plot = p12, filename = "./analysis/astrocytes/plots_final/subcortical_astros_regions_highlight.png", width = 10, height = 20, dpi = 600)

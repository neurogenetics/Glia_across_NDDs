library(scCustomize)
library(Seurat)
library(tidyverse)
library(svglite)
library(cowplot)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

micro <- readRDS("./analysis/microglia/seurat_objects/microglia_annotated.rds")

source("./code_organized/functions/plotting_colors.R")

##################################################
##################################################
##################################################

cluster_order <- rev(c("Micro_Homeo", 
                       "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                       "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                       "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN", 
                       "BAM", "Monocyte", "Lymphocyte"))

Idents(micro) <- factor(Idents(micro), levels = rev(cluster_order))

p1 <- DimPlot_scCustom(micro, repel = T, label.box = T, raster = F, colors_use = microglia_colors) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_fixed()

ggsave(p1, filename = "./analysis/microglia/plots_final/microglia_UMAP_clusters.png", height = 6.5, width = 9.5, dpi = 600)

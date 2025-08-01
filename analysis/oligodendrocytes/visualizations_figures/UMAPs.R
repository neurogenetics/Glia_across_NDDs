library(scCustomize)
library(Seurat)
library(tidyverse)
library(svglite)
library(cowplot)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

oligos <- readRDS("./analysis/oligodendrocytes/seurat_objects/oligos_clustered_ANNOTATED.rds")

source("./code_organized/functions/plotting_colors.R")

##################################################
##################################################
##################################################

cluster_order <- rev(c("Oligo_OPALIN", "Oligo_RBFOX1", "Oligo_CSMD1", "Oligo_CNTN1", "Oligo_Stress_HSPH1"))

Idents(oligos) <- factor(Idents(oligos), levels = rev(cluster_order))

p1 <- DimPlot_scCustom(oligos, repel = T, label.box = T, raster = F, colors_use = oligo_colors) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_fixed()

ggsave(p1, filename = "./analysis/oligodendrocytes/plots_final/oligos_UMAP_clusters.png", height = 6.5, width = 9.5, dpi = 600)

library(Seurat)
library(scCustomize)
library(tidyverse)
library(ggplot2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

micro <- readRDS("./analysis/microglia/seurat_objects/microglia_annotated.rds")
micro_signatures <- read.csv("./analysis/microglia/cluster_characterization/FGSEA/micro_signatures_to_test.csv", header = T)
micro_markers <- readRDS("./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea_ANNOTATED.rds")

##################################################

signatures_plot <- c("Sun_2023_MG0_Homeo", "Marshe_2025_CX3CR1.hi",
                     "Sun_2023_MG8_Inflamm2", "HUMICA_2025_c3_Ribo_DAM1",
                     "Sun_2023_MG4_Lipid", "Mancuso_2024_Xeno_DAM", "Dolan_2023_iMGL_c8_DAM",
                     "HUMICA_2025_c6_Lipo_DAM", "Marshe_2025_GPNMB.hi", 
                     "Silvin_2022_Human_DAMs", "Marshe_2025_APOE.hi", 
                     "Sun_2023_MG6_Stress", "Sun_2023_MG7_Glyco", "HUMICA_2025_c2_DIMs", "Marshe_2025_Stress",
                     "Chen_2024_PCDH9",
                     "Sun_2023_MG10_Inflamm3", "Marshe_2025_Chemokine",
                     "Sun_2023_MG5_Phago",
                     "Dolan_2023_iMGL_c10_Prolif", "Sun_2023_MG12_Cycling",
                     "Sun_2023_MG11_Antiviral", "Dolan_2023_iMGL_c11_Interferon")

micro_signatures <- micro_signatures[, colnames(micro_signatures) %in% signatures_plot]


micro_markers_filtered <- micro_markers[micro_markers$cluster == "Micro_DAM_Int2" &
                                          micro_markers$avg_log2FC > 0.5 &
                                          micro_markers$p_val_adj <= 0.05, ]


g1 <- intersect(micro_signatures$Sun_2023_MG11_Antiviral[1:100], micro_markers_filtered$gene[1:20])
g2 <- intersect(micro_signatures$Dolan_2023_iMGL_c11_Interferon[1:100], micro_markers_filtered$gene[1:20])
g3 <- intersect(micro_signatures$HUMICA_2025_c2_DIMs[1:100], micro_markers_filtered$gene[1:20])
g4 <- intersect(micro_signatures$Marshe_2025_Stress[1:100], micro_markers_filtered$gene[1:20])


gm <- g1
gm <- union(g1, g2)
gm <- union(g1, union(g2, g3))
gm <- union(g1, union(g2, union(g3, g4)))

##################################################

Idents(micro) <- "cluster_anno"

cluster_order <- rev(c("Micro_Homeo", 
                   "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                   "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                   "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN", 
                   "BAM", "Monocyte", "Lymphocyte"))

Idents(micro) <- factor(Idents(micro), levels = cluster_order)

genes_plot <- unique(c("CSF1R", "TMEM119", "AIF1", "CX3CR1", "P2RY12", "RASGEF1C", "SYNDIG1", 
                       "FOXP1", "APOE", "LRRK2", "DTNA",
                       "SPP1", "MYO1E", "PPARG", "GLDN",
                       "STARD13", "NHSL1", "CPM", "MITF",
                       "GPNMB", "CTSD", "CD9", "IQGAP2", "ABCA1", "PLIN2", 
                       "HSPA1A", "UBC", "P4HA1", "SLC2A3", "HIF1A", "DUSP1", "RGS1",
                       "PCDH9", "PLP1", "IL1RAPL1", "SLC44A1", 
                       "CD83", "CCL3", "CCL4", "NFKB1",
                       "CD163", "F13A1", 
                       "HELLS", "CENPP", "CLSPN", 
                       "IFI44L", "MX1", "MX2", "IFIT3",
                       "MRC1", "LYVE1",
                       "VCAN", "EMILIN2", "RIPOR2", 
                       "CD44", "CD72", "SDC2"))

p <- DotPlot_scCustom(micro, features = genes_plot, colors_use = c("#bbe7f2", "white", "#69078a")) + 
  RotatedAxis() + 
  theme(text = element_text(family = "Arial"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(p, filename = "./analysis/microglia/plots/dotplot_markers.png", width = 20, height = 6.2, dpi = 600)

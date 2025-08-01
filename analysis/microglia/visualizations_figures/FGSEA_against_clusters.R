library(tidyverse)
library(ggplot2)
library(scales)
library(svglite)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

fgsea_df <- readRDS("./analysis/microglia/cluster_characterization/FGSEA/micro_FGSEA_signatures_output_ANNOTATED.rds")

##################################################

micro_cluster_order <- rev(c("Micro_Homeo",  "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                             "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                             "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN"))

fgsea_df$cluster <- factor(fgsea_df$cluster, levels = micro_cluster_order)

##################################################

fgsea_df$signature <- gsub(".hi", "-hi", fgsea_df$signature)
fgsea_df$signature <- gsub("Keren.Shaul", "Keren-Shaul", fgsea_df$signature)
fgsea_df$signature <- gsub("APP.PS1", "APP-PS1", fgsea_df$signature)
fgsea_df$signature <- gsub("Martins.Ferreira", "Martins-Ferreira", fgsea_df$signature)

signatures_order <- c("Sun_2023_c0_Homeo", "Marshe_2025_CX3CR1-hi", "Marshe_2025_GRID2-hi",
                      "Martins-Ferreira_2025_c3_Ribo.DAM1", "Martins-Ferreira_2025_c5_Ribo.DAM2", "Mancuso_2024_Xenograft_Ribo",
                      "Gerrits_2021_c7_AD1", "Gerrits_2021_c9_AD1", "Gerrits_2021_c10_AD1",
                      "Sun_2023_c4_Lipid", "Martins-Ferreira_2025_c6_Lipo.DAM", "Marshe_2025_GPNMB-hi", "Silvin_2022_Human_DAMs",
                      "Dolan_2023_iMGL_c8_DAM", "Mancuso_2024_Xenograft_DAM", "Keren-Shaul_2017_5xFAD_DAM", "Krasemann_2017_APP-PS1_MGnD",
                      "Sun_2023_c6_Stress", "Sun_2023_c7_Glyco", "Marshe_2025_Senescence", "Marshe_2025_Stress", "Martins-Ferreira_2025_c2_DIMs",
                      "Mancuso_2024_Xenograft_CRM1", "Mancuso_2024_Xenograft_CRM2",
                      "Chen_2024_PCDH9-hi",
                      "Sun_2023_c10_Inflamm3",
                      "Sun_2023_c5_Phago",
                      "Sun_2023_c12_Cycling", "Dolan_2023_iMGL_c6_Prolif", "Dolan_2023_iMGL_c9_Prolif", "Dolan_2023_iMGL_c10_Prolif")

fgsea_df$signature <- factor(fgsea_df$signature, levels = signatures_order)

##################################################

fgsea_df$significant[is.na(fgsea_df$significant)] <- FALSE
fgsea_df <- na.omit(fgsea_df)

##################################################

fgsea_df$log10padj_squished <- pmin(pmax(fgsea_df$log10padj, 0), 50)
fgsea_df$NES_squished <- pmin(pmax(fgsea_df$NES, -3), 3)

##################################################

p1 <- ggplot(fgsea_df, aes(x = cluster, y = signature)) +
  geom_point(aes(size = log10padj_squished, fill = NES_squished, color = significant), shape = 21, stroke = 0.5) +
  scale_size(range = c(1, 10), limits = c(0, 50)) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") + 
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        text = element_text(family = "Arial"),
        axis.text.y = element_text(size = 12)) +
  labs(x = NULL, y = NULL, size = "-log10(padj)", color = "Significant", fill = "Normalized Enrichment") +
  coord_flip()

svglite("./analysis/microglia/plots_final/fgsea_against_clusters.svg", width = 14, height = 5.5)
plot(p1)
dev.off()

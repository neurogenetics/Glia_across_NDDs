library(tidyverse)
library(ggplot2)
library(svglite)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# CORTICAL ASTROS

##################################################

fgsea_df <- readRDS("./analysis/astrocytes/cluster_characterization/FGSEA/astro_FGSEA_signatures_output_ANNOTATED.rds")

##################################################

astro_cluster_order <- rev(c("Astro_Protoplasmic_GRM3",  "Astro_Fibrous_CD44", "Astro_Fibrous_GRIA1", 
                             "Astro_Reactive_SERPINA3", "Astro_Stress_HSPH1", "Astro_Stress_TXNRD1"))

fgsea_df$cluster <- factor(fgsea_df$cluster, levels = astro_cluster_order)

##################################################

fgsea_df$signature <- gsub("Serrano.Pozo", "Serrano-Pozo", fgsea_df$signature)

signatures_order <- c("Green_2024_c1_Homeostatic", "Green_2024_c2_Homeostatic", "Green_2024_c3_Mitophagy", "Green_2024_c4_Reactive", 
                      "Green_2024_c5_Reactive", "Green_2024_c6", "Green_2024_c7_Interferon", "Green_2024_c8_Stress", 
                      "Green_2024_c9_Stress", "Green_2024_c10_Stress", unique(fgsea_df$signature)[11:26])

fgsea_df$signature <- factor(fgsea_df$signature, levels = signatures_order)

##################################################

fgsea_df$significant[is.na(fgsea_df$significant)] <- FALSE
fgsea_df <- na.omit(fgsea_df)

##################################################

fgsea_df$log10padj_squished <- pmin(pmax(fgsea_df$log10padj, 0), 100)
fgsea_df$NES_squished <- pmin(pmax(fgsea_df$NES, -3), 3)

##################################################

p1 <- ggplot(fgsea_df, aes(x = cluster, y = signature)) +
  geom_point(aes(size = log10padj_squished, fill = NES_squished, color = significant), shape = 21, stroke = 0.5) +
  scale_size(range = c(1, 10), limits = c(0, 100)) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") + 
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        text = element_text(family = "Arial"),
        axis.text.y = element_text(size = 12)) +
  labs(x = NULL, y = NULL, size = "-log10(padj)", color = "Significant", fill = "Normalized Enrichment") +
  coord_flip()

svglite("./analysis/astrocytes/plots_final/fgsea_against_clusters.svg", width = 14, height = 3.8)
plot(p1)
dev.off()

##################################################
##################################################
##################################################

# ALL ASTROS

##################################################

fgsea_df <- readRDS("./analysis/astrocytes/cluster_characterization/FGSEA/ALL_astros_FGSEA_signatures_output_ANNOTATED.rds")

##################################################

astro_cluster_order <- rev(c("Astro_Protoplasmic_GRM3",  "Astro_Fibrous_DLCK1", "Astro_Fibrous_GRIA1", 
                             "Astro_KCND2", "Astro_Reactive_SERPINA3", "Astro_LAMA2", "Astro_Stress_HSPH1"))

fgsea_df$cluster <- factor(fgsea_df$cluster, levels = astro_cluster_order)

##################################################

fgsea_df$signature <- gsub("Serrano.Pozo", "Serrano-Pozo", fgsea_df$signature)

signatures_order <- c("Green_2024_c1_Homeostatic", "Green_2024_c2_Homeostatic", "Green_2024_c3_Mitophagy", "Green_2024_c4_Reactive", 
                      "Green_2024_c5_Reactive", "Green_2024_c6", "Green_2024_c7_Interferon", "Green_2024_c8_Stress", 
                      "Green_2024_c9_Stress", "Green_2024_c10_Stress", unique(fgsea_df$signature)[11:26])

fgsea_df$signature <- factor(fgsea_df$signature, levels = signatures_order)

##################################################

fgsea_df$significant[is.na(fgsea_df$significant)] <- FALSE
fgsea_df <- na.omit(fgsea_df)

##################################################

fgsea_df$log10padj_squished <- pmin(pmax(fgsea_df$log10padj, 0), 100)
fgsea_df$NES_squished <- pmin(pmax(fgsea_df$NES, -3), 3)

##################################################

p2 <- ggplot(fgsea_df, aes(x = cluster, y = signature)) +
  geom_point(aes(size = log10padj_squished, fill = NES_squished, color = significant), shape = 21, stroke = 0.5) +
  scale_size(range = c(1, 10), limits = c(0, 100)) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") + 
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        text = element_text(family = "Arial"),
        axis.text.y = element_text(size = 12)) +
  labs(x = NULL, y = NULL, size = "-log10(padj)", color = "Significant", fill = "Normalized Enrichment") +
  coord_flip()

svglite("./analysis/astrocytes/plots_final/ALL_astros_fgsea_against_clusters.svg", width = 14, height = 4.5)
plot(p2)
dev.off()

library(tidyverse)
library(ggplot2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

fgsea_df <- readRDS("./analysis/oligodendrocytes/cluster_characterization/FGSEA/oligo_FGSEA_signatures_output_ANNOTATED.rds")

##################################################

oligo_cluster_order <- rev(c("Oligo_OPALIN",  "Oligo_RBFOX1", "Oligo_CSMD1", 
                             "Oligo_CNTN1", "Oligo_Stress_HSPH1"))

fgsea_df$cluster <- factor(fgsea_df$cluster, levels = oligo_cluster_order)

##################################################

signatures_order <- c("Green_2024_c1", "Green_2024_c2", "Green_2024_c3", "Green_2024_c4", 
                      "Green_2024_c5", "Green_2024_c6", "Green_2024_c7", "Green_2024_c8", 
                      "Green_2024_c9", "Green_2024_c10", "Green_2024_c11", "Green_2024_c12",
                      unique(fgsea_df$signature)[13:18])

fgsea_df$signature <- factor(fgsea_df$signature, levels = signatures_order)

##################################################

fgsea_df$significant[is.na(fgsea_df$significant)] <- FALSE
fgsea_df <- na.omit(fgsea_df)

##################################################

fgsea_df$log10padj_squished <- pmin(pmax(fgsea_df$log10padj, 0), 50)
fgsea_df$NES_squished <- pmin(pmax(fgsea_df$NES, -2), 2)

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

svglite("./analysis/oligodendrocytes/plots_final/fgsea_against_clusters.svg", width = 9.5, height = 2.75)
plot(p1)
dev.off()

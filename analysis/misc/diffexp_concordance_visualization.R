library(tidyverse)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(svglite)

setwd("/data/ADRD/glia_across_NDDs")

########################################
########################################
########################################

# plots of concordance of nebula and DESeq results

astro_res_list <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_diffexp_results_for_mashr.rds")
micro_res_list <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_results_for_mashr.rds")
oligo_res_list <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_results_for_mashr.rds")


astro_spearman <- readRDS("./analysis/astrocytes/differential_expression/DESeq_nebula_correlation.rds")
micro_spearman <- readRDS("./analysis/microglia/differential_expression/DESeq_nebula_correlation.rds")
oligo_spearman <- readRDS("./analysis/oligodendrocytes/differential_expression/DESeq_nebula_correlation.rds")

source("./code_organized/functions/plotting_colors.R")

########################################
########################################
########################################

# stacked histogram of correlation values

spearman <- rbind(astro_spearman, micro_spearman, oligo_spearman)
spearman$dataset <- ifelse(grepl("AMP-PD", spearman$comparison), "NM_2024",
                           ifelse(grepl("Pineda", spearman$comparison), "Pineda_2024",
                                  ifelse(grepl("Mathys", spearman$comparison), "Mathys_2024", "Gerrits_2022")))


p1 <- ggplot(spearman, aes(x = spearman_zscore, fill = dataset)) +
  geom_histogram(binwidth = 0.05, boundary = -0.2) +
  labs(title = "Correlation between DESeq2 and nebula z-scored log2FC", 
       y = "Proportion",
       x = "Spearman coefficient") +
  scale_fill_manual(values = study_colors) + 
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

svglite(filename = "./analysis/microglia/plots_final/DESeq_nebula_comparison_histogram.svg", width = 7.3, height = 3.5)
plot(p1)
dev.off()

########################################
########################################
########################################

# representative plots

p2 <- ggplot(micro_res_list[["Mathys_TC_MTG_CR"]], aes(x = DESeq_zscore, y = nebula_zscore)) + 
  geom_hex(bins = 75) +
  stat_cor(aes(label = after_stat(r.label)),
           method = "spearman", 
           cor.coef.name = "rho",
           label.x.npc = "left",
           label.y.npc = "top", 
           size = 5) +
  labs(x = "DESeq2 z-score", 
       y = "nebula z-score",
       title = "Microglia, Mathys 2024, MTG, CR vs. HC") +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

svglite(filename = "./analysis/microglia/plots_final/DESeq_nebula_comparison_microglia.svg", width = 5.2, height = 4.2)
plot(p2)
dev.off()

########################################

p3 <- ggplot(astro_res_list[["Gerrits_FC_PFC_FTDGRN"]], aes(x = DESeq_zscore, y = nebula_zscore)) + 
  geom_hex(bins = 75) +
  stat_cor(aes(label = after_stat(r.label)),
           method = "spearman", 
           cor.coef.name = "rho",
           label.x.npc = "left",
           label.y.npc = "top", 
           size = 5) +
  labs(x = "DESeq2 z-score", 
       y = "nebula z-score",
       title = "Astrocytes, Gerrits 2022, FC, FTD-GRN vs. HC") +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

svglite(filename = "./analysis/microglia/plots_final/DESeq_nebula_comparison_astrocytes.svg", width = 5.2, height = 4.2)
plot(p3)
dev.off()

########################################

p4 <- ggplot(oligo_res_list[["Pineda_PFC_C9ALS"]], aes(x = DESeq_zscore, y = nebula_zscore)) + 
  geom_hex(bins = 75) +
  stat_cor(aes(label = after_stat(r.label)),
           method = "spearman", 
           cor.coef.name = "rho",
           label.x.npc = "left",
           label.y.npc = "top", 
           size = 5) +
  labs(x = "DESeq2 z-score", 
       y = "nebula z-score",
       title = "Oligodendrocytes, Pineda 2024, PFC, C9ALS vs. HC") +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

svglite(filename = "./analysis/microglia/plots_final/DESeq_nebula_comparison_oligos.svg", width = 5.2, height = 4.2)
plot(p4)
dev.off()

library(tidyverse)
library(mashr)
library(ggplot2)
library(ggh4x)
library(svglite)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################
##################################################
##################################################

micro <- readRDS("./analysis/microglia/cluster_proportions/micro_props_mashr_obj_1v1.rds")
astro <- readRDS("./analysis/astrocytes/cluster_proportions/cortical_astros_props_mashr_obj_canonical.rds")
oligo <- readRDS("./analysis/oligodendrocytes/cluster_proportions/oligo_props_mashr_obj_canonical.rds")

source("./code_organized/functions/plotting_colors.R")

##################################################
##################################################
##################################################

micro_lfsr <- as.data.frame(micro$result$lfsr)
micro <- as.data.frame(micro$result$PosteriorMean)

micro_lfsr <- micro_lfsr %>%
  rownames_to_column(var = "celltype") %>%
  pivot_longer(cols = -all_of("celltype"), names_to = "d_r", values_to = "lfsr")

micro <- micro %>%
  rownames_to_column(var = "celltype") %>%
  pivot_longer(cols = -all_of("celltype"), names_to = "d_r", values_to = "beta") %>%
  left_join(micro_lfsr, by = c("celltype", "d_r"))

#########################

astro_lfsr <- as.data.frame(astro$result$lfsr)
astro <- as.data.frame(astro$result$PosteriorMean)

astro_lfsr <- astro_lfsr %>%
  rownames_to_column(var = "celltype") %>%
  pivot_longer(cols = -all_of("celltype"), names_to = "d_r", values_to = "lfsr")

astro <- astro %>%
  rownames_to_column(var = "celltype") %>%
  pivot_longer(cols = -all_of("celltype"), names_to = "d_r", values_to = "beta") %>%
  left_join(astro_lfsr, by = c("celltype", "d_r"))

#########################

oligo_lfsr <- as.data.frame(oligo$result$lfsr)
oligo <- as.data.frame(oligo$result$PosteriorMean)

oligo_lfsr <- oligo_lfsr %>%
  rownames_to_column(var = "celltype") %>%
  pivot_longer(cols = -all_of("celltype"), names_to = "d_r", values_to = "lfsr")

oligo <- oligo %>%
  rownames_to_column(var = "celltype") %>%
  pivot_longer(cols = -all_of("celltype"), names_to = "d_r", values_to = "beta") %>%
  left_join(oligo_lfsr, by = c("celltype", "d_r"))

#########################

micro <- separate(micro, col = "d_r", into = c("disease", "region"), sep = "_")
micro$ct = "Microglia"

astro <- separate(astro, col = "d_r", into = c("disease", "region"), sep = "_")
astro$ct = "Astrocytes"

oligo <- separate(oligo, col = "d_r", into = c("disease", "region"), sep = "_")
oligo$ct = "Oligodendrocytes"


t <- rbind(micro, astro, oligo)

##################################################

t$sig <- ifelse(t$beta > 0 & t$lfsr < 0.05, "enriched, lfsr < 0.05",
                ifelse(t$beta < 0 & t$lfsr < 0.05, "depleted, lfsr < 0.05", "n.s."))

##################################################

t$disease <- factor(t$disease, levels = c("ALS", "C9ALS", "FTLD", "C9FTLD",
                                          "FTD-GRN", "AD", "CR", "PD"))
disease_levels <- levels(factor(t$disease))
strip_colors <- disease_colors[disease_levels]


t$ct <- factor(t$ct, levels = c("Oligodendrocytes", "Microglia", "Astrocytes"))

##################################################

p1 <- ggplot(t, aes(x = ct, y = beta)) +
  geom_boxplot(width = 0.6, outliers = F) + 
  geom_jitter(width = 0.1, aes(colour = sig), alpha = 0.7) +
  scale_colour_manual(values = c("#05409e", "#e64a02", "grey40")) + 
  theme_bw() +
  labs(x = NULL,
       y = "βmashr") +
  ggh4x::facet_wrap2(~ disease, nrow = 2,
                     strip = strip_themed(background_x = elem_list_rect(fill = strip_colors))) +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(face = "bold", size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
  coord_flip()

svglite("./analysis/misc/props_summary_fig.svg", width = 9.5, height = 3.2)
plot(p1)
dev.off()

library(corrplot)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggh4x)

setwd("/data/ADRD/glia_across_NDDs")

########################################
########################################
########################################

source("./code_organized/functions/plot_proportions.R")
source("./code_organized/functions/plotting_colors.R")
source("./code_organized/functions/rename_datasets.R")

########################################
########################################
########################################

# CORTICAL

########################################

props_mashr <- readRDS("./analysis/astrocytes/cluster_proportions/cortical_astros_props_mashr_obj_canonical.rds")
props_raw <- readRDS("./analysis/astrocytes/cluster_proportions/cortical_astro_raw_donor_region_props.rds")

########################################

# heatmap of beta values with lfsr significance stars

col_fun2 = colorRamp2(c(-1, 0, 3), c("#05409e", "white", "#e64a02"))

fdr_fun1 <- function(j, i, x, y, w, h, fill) {
  if(props_mashr$result$lfsr[i, j] <= 0.0001) {
    grid.text("﹡﹡﹡﹡", x, y)
  } else if(props_mashr$result$lfsr[i, j] <= 0.001) {
    grid.text("﹡﹡﹡", x, y)
  } else if(props_mashr$result$lfsr[i, j] <= 0.01) {
    grid.text("﹡﹡", x, y)
  } else if(props_mashr$result$lfsr[i, j] <= 0.05) {
    grid.text("﹡", x, y)
  }
}

########################################

# make heatmap annotation

annotation_data <- data.frame(r_d <- colnames(props_mashr$result$PosteriorMean))
annotation_data$disease <- sub("_.*", "", annotation_data$r_d....colnames.props_mashr.result.PosteriorMean.)
annotation_data$region <- sub("^[^_]*_", "", annotation_data$r_d....colnames.props_mashr.result.PosteriorMean.)
annotation_data$dataset <- ifelse(grepl("ALS|FTLD", annotation_data$disease), "Pineda_2024",
                                  ifelse(grepl("AD|CR", annotation_data$disease), "Mathys_2024",
                                         ifelse(grepl("PD", annotation_data$disease), "NM_2024", "Gerrits_2022")))
annotation_data <- annotation_data %>%
  column_to_rownames(var = "r_d....colnames.props_mashr.result.PosteriorMean.")


ha = HeatmapAnnotation(Disease = annotation_data$disease,
                       `Brain region` = annotation_data$region, 
                       Dataset = annotation_data$dataset,
                       col = list(
                         Disease = disease_colors,
                         `Brain region` = brain_region_colors,
                         Dataset = study_colors),
                       gp = gpar(col = "black"),
                       annotation_legend_param = list(direction = "horizontal", title_position = "topleft"))

####################

postmean <- props_mashr$result$PosteriorMean

colnames(postmean) <- rename_datasets(colnames(postmean))

ht = Heatmap(as.matrix(postmean),
             col = col_fun2,
             cell_fun = fdr_fun1,
             rect_gp = gpar(col = "black"),
             name = "β_mashr",
             bottom_annotation = ha)

png('./analysis/astrocytes/plots_final/cortical_props_heatmap.png', res = 600, width = 14.5, height = 5.5, units = "in")
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = T)
dev.off()

########################################
########################################
########################################

# boxplots of raw props values

lfsr_df <- as.data.frame(props_mashr$result$lfsr)

####################

plot_props <- function(props_df, lfsr_df = NULL, dataset, clusters, region, disease, title){
  props_plot <- props_df[props_df$study == dataset &
                           props_df$cluster %in% clusters &
                           props_df$region == region & 
                           props_df$group %in% c(disease, "HC"),]
  
  d_r <- paste0(disease, "_", region)
  lfsr_plot <- lfsr_df[c(d_r)] %>%
    rownames_to_column(var = "cluster")
  lfsr_plot <- lfsr_plot[lfsr_plot$cluster %in% clusters, ]
  
  props_plot <- left_join(props_plot, lfsr_plot, by = "cluster") %>%
    dplyr::rename(lfsr = d_r) %>%
    group_by(cluster) %>%
    mutate(max_cluster_proportion = max(cluster_proportion)) %>%
    ungroup()
  
  props_plot$group <- factor(props_plot$group, levels = c("HC", disease))
  props_plot$cluster <- factor(props_plot$cluster, levels = clusters)
  
  cluster_levels <- levels(factor(props_plot$cluster))
  strip_colors <- astrocyte_colors[cluster_levels]

  p <- ggplot(props_plot, aes(x = group, y = cluster_proportion, fill = group)) + 
    geom_boxplot() + 
    scale_fill_manual(values = disease_colors) + 
    ggh4x::facet_wrap2(~ cluster, nrow = 1, scales = "free_y",
                       strip = strip_themed(background_x = elem_list_rect(fill = strip_colors))) +
    geom_point() + 
    geom_text(aes(x = 1.5, y = (max_cluster_proportion + (0.05*max_cluster_proportion)), 
                  label = paste0("lfsr = ", formatC(lfsr, format = "e", digits = 2)))) + 
    theme_bw() +
    labs(y = "Cluster proportion",
         title = title) + 
    theme(text = element_text(family = "Arial"),
          legend.position = "none",
          axis.title.x = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  return(p)
}

####################

p1 <- plot_props(props_df = props_raw, lfsr_df = lfsr_df, dataset = "Mathys", 
                 clusters = c("Astro_Fibrous_CD44", "Astro_Fibrous_GRIA1"), 
                 region = "EC", disease = "AD", title = "Mathys 2024, EC, AD vs. HC")

p2 <- plot_props(props_df = props_raw, lfsr_df = lfsr_df, dataset = "Pineda", 
                 clusters = c("Astro_Fibrous_CD44", "Astro_Fibrous_GRIA1"), 
                 region = "FC/PFC", disease = "C9FTLD", title = "Pineda 2024, PFC, C9FTLD vs. HC") +
  theme(axis.title.y = element_blank())

p3 <- plot_props(props_df = props_raw, lfsr_df = lfsr_df, dataset = "Gerrits", 
                 clusters = c("Astro_Fibrous_CD44", "Astro_Fibrous_GRIA1"), 
                 region = "FC/PFC", disease = "FTD-GRN", title = "Gerrits 2022, FC, FTD-GRN vs. HC") +
  theme(axis.title.y = element_blank())


p4 <- plot_grid(p1, p2, p3, nrow = 1)

svglite("./analysis/astrocytes/plots_final/raw_props_boxplots.svg", width = 15, height = 3.5)
plot(p4)
dev.off()

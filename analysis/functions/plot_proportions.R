library(tidyverse)
library(ggplot2)
library(ggh4x)

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
  
  p <- ggplot(props_plot, aes(x = group, y = cluster_proportion, fill = group)) + 
    geom_boxplot() + 
    scale_fill_manual(values = disease_colors) + 
    ggh4x::facet_wrap2(~ cluster, nrow = 1, scales = "free_y",
                       strip = strip_themed(background_x = elem_list_rect(fill = microglia_colors))) +
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

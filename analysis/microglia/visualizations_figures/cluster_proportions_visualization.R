library(corrplot)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(svglite)
library(mashr)

setwd("/data/ADRD/glia_across_NDDs")

########################################

source("./code_organized/functions/plot_proportions.R")
source("./code_organized/functions/plotting_colors.R")
source("./code_organized/functions/rename_datasets.R")

props_mashr <- readRDS("./analysis/microglia/cluster_proportions/micro_props_mashr_obj_1v1.rds")
props_raw <- readRDS("./analysis/microglia/cluster_proportions/micro_raw_donor_region_props.rds")

########################################
########################################
########################################

# plot of contributions

pl = data.frame(mi = get_estimated_pi(props_mashr))
pl$group = rownames(pl)

pl$group <- rename_datasets(pl$group)

pl <- pl[!grepl("2024|2022", pl$group), ]
pl <- t(pl)
pl <- as.data.frame(pl)
pl$individual_effects_sum <- c(0, "indiv_sum")
pl <- as.data.frame(t(pl))

pl <- pl %>%
  arrange(desc(mi))

pl$group <- factor(pl$group, levels = pl$group)



p1 <- ggplot(pl, aes(x = group, y = as.double(mi), fill = group)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c("ED_tPCA" = "red")) + 
  theme_bw() +
  labs(y = "Mixture proportion") + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(angle = 45, size = 12, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none")

svglite("./analysis/microglia/plots_final/proportions_mashr_mixtures_summary.svg", width = 6.6, height = 2.7)
plot(p1)
dev.off()

########################################
########################################
########################################

# mashr tPCA correlation plot

tPCA <- props_mashr$fitted_g$Ulist$ED_tPCA

colnames(tPCA) <- c("Pineda_2024_PFC_ALS", "Pineda_2024_PFC_C9ALS", "Pineda_2024_PFC_FTLD", "Pineda_2024_PFC_C9FTLD", 
                    "Pineda_2024_M1_ALS", "Pineda_2024_M1_C9ALS", "Pineda_2024_M1_FTLD", "Pineda_2024_M1_C9FTLD", 
                    "Mathys_2024_AnG_AD", "Mathys_2024_AnG_CR", "Mathys_2024_EC_AD", "Mathys_2024_EC_CR", 
                    "Mathys_2024_PFC_AD", "Mathys_2024_PFC_CR", "Mathys_2024_HIP_AD", "Mathys_2024_HIP_CR", 
                    "Mathys_2024_MTG_AD", "Mathys_2024_MTG_CR", "Mathys_2024_TH_AD", "Mathys_2024_TH_CR", 
                    "Gerrits_2022_FC_FTD-GRN", "Gerrits_2022_OC_FTD-GRN", "Gerrits_2022_TC_FTD-GRN", 
                    "NM_2024_DMV_PD", "NM_2024_PFC_PD", "NM_2024_GPi_PD", "NM_2024_M1_PD", "NM_2024_V1_PD")

rownames(tPCA) <- colnames(tPCA)

# bubble plot
col <- colorRampPalette(c("#05409e", "white", "#e64a02"))(200)

png('./analysis/microglia/plots_final/props_mashr_corrplot.png', res = 600, width = 10, height = 8, units = "in")
corrplot(tPCA, order = "hclust", col = col, tl.col = "black")
dev.off()

########################################
########################################
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

####################

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

colnames(postmean) <- colnames(tPCA)

ht = Heatmap(as.matrix(postmean),
             col = col_fun2,
             cell_fun = fdr_fun1,
             rect_gp = gpar(col = "black"),
             name = "β_mashr",
             bottom_annotation = ha)

png('./analysis/microglia/plots_final/props_heatmap.png', res = 600, width = 15, height = 6.5, units = "in")
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = T)
dev.off()

########################################

# boxplots of raw props values

lfsr_df <- as.data.frame(props_mashr$result$lfsr)

######

p1 <- plot_props(props_df = props_raw, lfsr_df = lfsr_df, dataset = "Gerrits", 
                 clusters = c("Micro_Homeo", "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB"), 
                 region = "FC/PFC", disease = "FTD-GRN", title = "Gerrits 2022, FC, FTD-GRN vs. HC")

p2 <- plot_props(props_df = props_raw, lfsr_df = lfsr_df, dataset = "Pineda", 
                 clusters = c("Micro_Homeo", "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB"), 
                 region = "M1", disease = "FTLD", title = "Pineda 2024, M1, FTLD vs. HC")

p3 <- plot_props(props_df = props_raw, lfsr_df = lfsr_df, dataset = "Mathys",
                 clusters = c("Micro_Homeo", "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB"),
                 region = "EC", disease = "AD", title = "Mathys 2024, EC, AD vs. HC")

# ggsave(plot = p1, filename = "./analysis/microglia/plots/Gerrits_FC_micro_props.png", width = 9.8, height = 3.3, dpi = 600)
# ggsave(plot = p2, filename = "./analysis/microglia/plots/Pineda_FC_micro_props.png", width = 9.8, height = 3.3, dpi = 600)

svglite(p1, file = "./analysis/microglia/plots_final/raw_props_Gerrits_FC.svg", width = 9.8, height = 3.3)
plot(p1)
dev.off()

svglite(p2, file = "./analysis/microglia/plots_final/raw_props_Pineda_M1_FTLD.svg", width = 9.8, height = 3.3)
plot(p2)
dev.off()

svglite(p3, file = "./analysis/microglia/plots_final/raw_props_Mathys_EC_AD.svg", width = 9.8, height = 3.3)
plot(p3)
dev.off()

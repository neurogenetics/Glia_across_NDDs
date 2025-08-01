library(corrplot)
library(ComplexHeatmap)
library(tidyverse)
library(mashr)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(enrichplot)
library(svglite)

setwd("/data/ADRD/glia_across_NDDs")

########################################

micro_mashr <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_mashr_object.rds")
tPCA_betas <- readRDS("./analysis/microglia/differential_expression/microglia_mashr_tPCA_betas.rds")
tPCA_lfsrs <- readRDS("./analysis/microglia/differential_expression/microglia_mashr_tPCA_lfsrs.rds")

gobp_c1_c2 <- readRDS("./analysis/microglia/differential_expression/microglia_gobp_c1_c2.rds")
gobp_c4_c6 <- readRDS("./analysis/microglia/differential_expression/microglia_gobp_c4_c6.rds")
gobp_all <- readRDS("./analysis/microglia/differential_expression/microglia_gobp_all.rds")

source("./code_organized/functions/plotting_colors.R")
source("./code_organized/functions/rename_datasets.R")

########################################
########################################
########################################

# plot of contributions

pl = data.frame(mi = get_estimated_pi(micro_mashr))
pl$group = rownames(pl)
pl <- pl %>%
  arrange(desc(mi))

pl$group <- rename_datasets(pl$group)

pl$group <- factor(pl$group, levels = pl$group)



p1 <- ggplot(pl, aes(x = group, y = mi)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  labs(y = "Mixture proportion") + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank())

ggsave(plot = p1, filename = "./analysis/microglia/plots_final/diffexp_mashr_mixtures.png", width = 9.5, height = 5.5, dpi = 600)

########################################
########################################
########################################

# corrplot of different covariance matrices

tPCA_corr <- micro_mashr$fitted_g$Ulist$ED_tPCA
colnames(tPCA_corr) <- rename_datasets(colnames(tPCA_corr))
rownames(tPCA_corr) <- rename_datasets(colnames(tPCA_corr))

col <- colorRampPalette(c("#05409e", "white", "#e64a02"))(200)

png('./analysis/microglia/plots_final/diffexp_mashr_corrplot.png', res = 600, width = 10, height = 8, units = "in")
corrplot(tPCA_corr, order = "hclust", col = col, tl.col = "black")
dev.off()

########################################
########################################
########################################

# make k-means clustered heatmap

colnames(tPCA_betas) <- rename_datasets(colnames(tPCA_betas))


annotation_data <- data.frame(r_d <- colnames(tPCA_betas))
annotation_data$disease <- str_split(annotation_data$r_d....colnames.tPCA_betas., "_", simplify = TRUE)[, 3]
annotation_data$region <- str_split(annotation_data$r_d....colnames.tPCA_betas., "_", simplify = TRUE)[, 4]
annotation_data$dataset <- ifelse(grepl("ALS|FTLD", annotation_data$disease), "Pineda_2024",
                                  ifelse(grepl("AD|CR", annotation_data$disease), "Mathys_2024",
                                         ifelse(grepl("PD", annotation_data$disease), "NM_2024", "Gerrits_2022")))
annotation_data <- annotation_data %>%
  column_to_rownames(var = "r_d....colnames.tPCA_betas.")

ha = HeatmapAnnotation(Disease = annotation_data$disease,
                       `Brain region` = annotation_data$region, 
                       Dataset = annotation_data$dataset,
                       col = list(
                         Disease = disease_colors,
                         `Brain region` = brain_region_colors,
                         Dataset = study_colors),
                       gp = gpar(col = "black"),
                       annotation_name_side = "left")

col_fun <- colorRamp2(c(-1, 0, 1), c("#05409e", "white", "#e64a02"))

set.seed(12345)
hmp1 <- Heatmap(as.matrix(tPCA_betas),
                show_row_names = F,
                cluster_rows = F,
                row_km = 6, 
                col = col_fun,
                name = "β_mashr",
                bottom_annotation = ha,
                use_raster = F)
hmp1 <- draw(hmp1, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = T)


png('./analysis/microglia/plots_final/diffexp_hmp.png', res = 600, width = 11.5, height = 6.5, units = "in")
draw(hmp1)
dev.off()

########################################
########################################
########################################

# plotting gene ontology stuff

gobp_c1_c2 <- readRDS("./analysis/microglia/differential_expression/microglia_gobp_c1_c2.rds")
gobp_c4_c6 <- readRDS("./analysis/microglia/differential_expression/microglia_gobp_c4_c6.rds")

########################################

# cluster 1/2
gobp_plot <- gobp_c1_c2[1:5,]
gobp_plot$`-log10(FDR)` <- -log10(gobp_plot$p.adjust)
gobp_plot <- gobp_plot %>%
  arrange(desc(`-log10(FDR)`))
gobp_plot$Description <- factor(gobp_plot$Description, levels = rev(gobp_plot$Description))

p1 <- ggplot(gobp_plot, aes(x = `-log10(FDR)`, y = Description)) +
  geom_col(fill = "#3171c4") +
  theme_bw() + 
  labs(title = "K-means cluster 1, 2",
       x = "-log10(FDR)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none",
        text = element_text(family = "Arial")) +
  geom_vline(xintercept = (-log10(0.05)), linetype = "dashed", color = "red", size = 2)



# cluster 4/6
gobp_plot <- gobp_c4_c6[1:5,]
gobp_plot$`-log10(FDR)` <- -log10(gobp_plot$p.adjust)
gobp_plot <- gobp_plot %>%
  arrange(desc(`-log10(FDR)`))
gobp_plot$Description <- c(gobp_plot$Description[1:3], "cellular response to peptide\nhormone stimulus", gobp_plot$Description[5])
gobp_plot$Description <- factor(gobp_plot$Description, levels = rev(gobp_plot$Description))

p2 <- ggplot(gobp_plot, aes(x = `-log10(FDR)`, y = Description)) +
  geom_col(fill = "#3171c4") +
  theme_bw() + 
  labs(title = "K-means clusters 4, 6",
       x = "-log10(FDR)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none",
        text = element_text(family = "Arial")) +
  geom_vline(xintercept = (-log10(0.05)), linetype = "dashed", color = "red", size = 2)

####################

svglite(filename = "./analysis/microglia/plots_final/GOBP_c1_c2.svg", width = 7, height = 3.2)
plot(p1)
dev.off()

svglite(filename = "./analysis/microglia/plots_final/GOBP_c4_c6.svg", width = 7, height = 3.2)
plot(p2)
dev.off()

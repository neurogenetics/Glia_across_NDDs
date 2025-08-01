library(ComplexHeatmap)
library(tidyverse)
library(circlize)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

DAM_TFs <- readRDS("./analysis/microglia/DAM_signature_discovery/tf_varpart/final_summary_DAM_TFs.rds")
other_TFs <- readRDS("./analysis/microglia/DAM_signature_discovery/tf_varpart/other_TF_contributions_summary.rds")

##################################################
##################################################
##################################################

# organize data & FDR-correct

other_TFs <- other_TFs[c("tf_var_expl", "TF", "test", "anova")]
colnames(other_TFs) <- colnames(DAM_TFs)

df_merged <- rbind(DAM_TFs, other_TFs)

df_merged$padj <- p.adjust(df_merged$pval, method = "BH")



anova <- df_merged %>%
  dplyr::select(TF, dataset_region, padj) %>%
  pivot_wider(names_from = TF, values_from = padj) %>%
  column_to_rownames(var = "dataset_region") %>%
  t(.) %>%
  as.matrix(.)

colnames(anova) <- c("Mathys_2024_HIP", "Mathys_2024_EC", "Mathys_2024_TH", "Mathys_2024_PFC", "Mathys_2024_AnG", "Mathys_2024_MTG", 
                     "NM_2024_DMV", "NM_2024_GPi", "NM_2024_M1", "NM_2024_PFC", "NM_2024_V1",
                     "Gerrits_2022_FC", "Gerrits_2022_TC", "Gerrits_2022_OC", "Pineda_2024_M1", "Pineda_2024_PFC")
anova <- anova[, c("Mathys_2024_TH", "Mathys_2024_EC", "Mathys_2024_MTG", "Mathys_2024_HIP", "Mathys_2024_AnG", "Mathys_2024_PFC", 
                   "NM_2024_DMV", "NM_2024_GPi", "NM_2024_M1", "NM_2024_PFC", "NM_2024_V1",
                   "Gerrits_2022_TC", "Gerrits_2022_OC", "Gerrits_2022_FC", "Pineda_2024_M1", "Pineda_2024_PFC")]
anova <- anova[c("MITF", "PPARG", "ARID5B", "MAFB", "combined", "CEBPA", "SPI1"),]
anova[is.na(anova)] <- 1

saveRDS(anova, file = "./analysis/microglia/DAM_signature_discovery/tf_varpart/anova_res_ALLcomps_df.rds")

##################################################

variance <- df_merged %>%
  dplyr::select(TF, dataset_region, var_explained) %>%
  pivot_wider(names_from = TF, values_from = var_explained) %>%
  column_to_rownames(var = "dataset_region") %>%
  t(.) %>%
  as.matrix(.)

colnames(variance) <- c("Mathys_2024_HIP", "Mathys_2024_EC", "Mathys_2024_TH", "Mathys_2024_PFC", "Mathys_2024_AnG", "Mathys_2024_MTG", 
                        "NM_2024_DMV", "NM_2024_GPi", "NM_2024_M1", "NM_2024_PFC", "NM_2024_V1",
                        "Gerrits_2022_FC", "Gerrits_2022_TC", "Gerrits_2022_OC", "Pineda_2024_M1", "Pineda_2024_PFC")
variance <- variance[, colnames(anova)]
variance <- variance[c("MITF", "PPARG", "ARID5B", "MAFB", "combined", "CEBPA", "SPI1"),]
rownames(variance) <- c("MITF", "PPARG", "ARID5B", "MAFB", "All 4\n(combinatorial)", "CEBPA", "SPI1")

saveRDS(variance, file = "./analysis/microglia/DAM_signature_discovery/tf_varpart/variance_explained_ALLcomps_df.rds")

##################################################
##################################################
##################################################

col_fun = colorRamp2(breaks = c(0, 1), colors = c("white", "red"))

##################################################

fdr_fun <- function(j, i, x, y, w, h, fill) {
  if(anova[i, j] <= 0.0001) {
    grid.text("﹡﹡﹡﹡", x, y)
  } else if(anova[i, j] <= 0.001) {
    grid.text("﹡﹡﹡", x, y)
  } else if(anova[i, j] <= 0.01) {
    grid.text("﹡﹡", x, y)
  } else if(anova[i, j] <= 0.05) {
    grid.text("﹡", x, y)
  }
}

row_groups <- c(rep("Group1", 4), rep("Group2", 1), rep("Group3", 2))


ht = Heatmap(variance,
             col = col_fun,
             cluster_rows = F,
             cluster_columns = F,
             cell_fun = fdr_fun,
             rect_gp = gpar(col = "black"),
             row_split = row_groups,
             name = "Proportion variance explained",
             heatmap_legend_param = list(title_position = "topleft",
                                         legend_direction = "horizontal",
                                         legend_width = unit(3, "cm"),
                                         title_gp = gpar(fontsize = 12)),
             row_title = NULL)

png("./analysis/microglia/plots_final/TF_varpart_hmp.png", width = 7, height = 4.5, units = "in", res = 600)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

library(ggplot2)
library(ggdendro)
library(tidyverse)
library(Seurat)
library(cowplot) 
library(extrafont)
library(scCustomize)
library(patchwork)
library(svglite)
library(scales)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################
##################################################
##################################################

source("./code_organized/functions/plotting_colors.R")

theme_empty <- function(x) {
  theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                          axis.text.y = element_blank(), axis.title = element_blank(), 
                          legend.position = "none")
}

##################################################
##################################################
##################################################

#############
# MICROGLIA #
#############

micro <- readRDS("./analysis/microglia/seurat_objects/microglia_annotated.rds")
micro_meta <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")

micro_cluster_order <- rev(c("Micro_Homeo",  "Micro_DAM_Int1", "Micro_DAM_SPP1", "Micro_DAM_Int2", "Micro_DAM_GPNMB",
                             "Micro_Inflamm_Stress", "Micro_Inflamm_PCDH9", "Micro_Inflamm_CD83",
                             "Micro_Phago_CD163", "Micro_Prolif", "Micro_IFN", 
                             "BAM", "Monocyte", "Lymphocyte"))

##################################################

# number of cells, plotted by dataset

micro_ncells_data <- micro_meta %>%
  group_by(cluster_anno, dataset) %>%
  summarize(ncells = n()) %>%
  ungroup()

micro_ncells_data$dataset <- ifelse(micro_ncells_data$dataset == "Gerrits", paste0(micro_ncells_data$dataset, "_2022"),
                                    ifelse(micro_ncells_data$dataset == "AMP-PD", "NM_2024", paste0(micro_ncells_data$dataset, "_2024")))

micro_ncells_plot <- ggplot(micro_ncells_data, aes(x = ncells, y = cluster_anno, fill = dataset)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = study_colors) +
  scale_y_discrete(limits = micro_cluster_order) +
  scale_x_continuous(position = "top", breaks = c(0, 50000, 100000, 150000), labels = label_comma()) +
  theme_empty() +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 10, angle = 60, hjust = 0), 
        axis.ticks.x = element_line(linewidth = 0.2), 
        axis.line.x = element_line(linewidth = 0.2),
        axis.line.y = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2), 
        text = element_text(family = "Arial"))

##################################################

# proportions of cells cortical/subcortical

micro_meta$broad_region <- ifelse(micro_meta$region %in% c("DMV", "GPi", "TH"), "Subcortical", "Cortical")

micro_region_data <- micro_meta %>%
  group_by(cluster_anno, broad_region) %>%
  summarise(n = n()) %>% 
  group_by(broad_region) %>%
  mutate(num_cells_in_region = sum(n)) %>%
  mutate(prop = n / num_cells_in_region * 100)

micro_region_data$broad_region <- factor(micro_region_data$broad_region, levels = c("Subcortical", "Cortical"))

micro_region_plot <- ggplot(micro_region_data, aes(x = prop, y = cluster_anno, fill = broad_region)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  scale_fill_manual(values = c(Cortical = "#9e1f63", Subcortical = "#6fc251")) +
  scale_y_discrete(limits = micro_cluster_order) +
  scale_x_continuous(position = "top", transform = "log1p") +
  theme_empty() +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 10, angle = 60, hjust = 0), 
        axis.ticks.x = element_line(linewidth = 0.2), axis.line.x = element_line(linewidth = 0.2),
        text = element_text(family = "Arial")) 

##################################################

# dotplot of marker genes

micro_genes <- unique(c("CSF1R", "TMEM119", "AIF1", "CX3CR1", "P2RY12", "RASGEF1C", "SYNDIG1", 
                        "FOXP1", "APOE", "LRRK2", "DTNA",
                        "SPP1", "MYO1E", "PPARG", "GLDN",
                        "STARD13", "NHSL1", "CPM", "MITF",
                        "GPNMB", "CTSD", "CD9", "IQGAP2", "ABCA1", "PLIN2", 
                        "HSPA1A", "UBC", "P4HA1", "SLC2A3", "HIF1A", "DUSP1", "RGS1",
                        "PCDH9", "PLP1", "IL1RAPL1", "SLC44A1", 
                        "CD83", "CCL3", "CCL4", "NFKB1",
                        "CD163", "F13A1", 
                        "HELLS", "CENPP", "CLSPN", 
                        "IFI44L", "MX1", "MX2", "IFIT3",
                        "MRC1", "LYVE1",
                        "VCAN", "EMILIN2", "RIPOR2", 
                        "CD44", "CD72", "SDC2"))

micro_dotplot <- DotPlot(micro, features = micro_genes, assay = "sketch", dot.min = 0.05, dot.scale = 6) +
  aes(colour = avg.exp.scaled) +
  scale_color_gradient(limits = c(-1, 2), low = "#d9f4fa", high = "#64288a", oob = squish) +
  scale_x_discrete(limits = micro_genes, position = "top") +
  scale_y_discrete(limits = micro_cluster_order) +
  theme(legend.position = "right", axis.title = element_blank(), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.2), axis.ticks = element_line(linewidth = 0.2), 
        axis.text.x = element_text(size = 10, angle = 60, hjust = 0), axis.text.y = element_blank(),
        text = element_text(family = "Arial")) +
  NoLegend()

########

# combine plots

p1 <- plot_grid(micro_ncells_plot, micro_region_plot, micro_dotplot, nrow = 1, ncol = 3, rel_widths = c(1.5, 1.5, 7),
                align = "h", axis = "tb") %>% print()

########

# print a version w/ legend

svglite(filename = "./analysis/misc/micro_dotplot_with_legend.svg", width = 15.5, height = 9)
plot(micro_dotplot)
dev.off()

##################################################
##################################################
##################################################

##############
# ASTROCYTES #
##############

astro <- readRDS("./analysis/astrocytes/seurat_objects/all_astrocytes_clustered_ANNOTATED.rds")
astro_meta <- readRDS("./analysis/astrocytes/all_astros_metadata_ANNOTATED.rds")

astro_cluster_order <- rev(c("Astro_Protoplasmic_GRM3", "Astro_Fibrous_DLCK1", "Astro_Fibrous_GRIA1",
                             "Astro_KCND2", "Astro_Reactive_SERPINA3", "Astro_LAMA2", "Astro_Stress_HSPH1"))

##################################################

# number of cells, plotted by dataset

astro_ncells_data <- astro_meta %>%
  group_by(cluster_anno, dataset) %>%
  summarize(ncells = n()) %>%
  ungroup()

astro_ncells_data$dataset <- ifelse(astro_ncells_data$dataset == "Gerrits", paste0(astro_ncells_data$dataset, "_2022"),
                                    ifelse(astro_ncells_data$dataset == "AMP-PD", "NM_2024", paste0(astro_ncells_data$dataset, "_2024")))

astro_ncells_plot <- ggplot(astro_ncells_data, aes(x = ncells, y = cluster_anno, fill = dataset)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = study_colors) +
  scale_y_discrete(limits = astro_cluster_order) +
  scale_x_continuous(position = "top", labels = label_comma()) +
  theme_empty() +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 10, angle = 60, hjust = 0), 
        axis.ticks.x = element_line(linewidth = 0.2), 
        axis.line.x = element_line(linewidth = 0.2),
        axis.line.y = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2), 
        text = element_text(family = "Arial"))

##################################################

# proportions of cells cortical/subcortical

astro_meta$broad_region <- ifelse(astro_meta$region %in% c("DMV", "GPi", "TH"), "Subcortical", "Cortical")

astro_region_data <- astro_meta %>%
  group_by(cluster_anno, broad_region) %>%
  summarise(n = n()) %>% 
  group_by(broad_region) %>%
  mutate(num_cells_in_region = sum(n)) %>%
  mutate(prop = n / num_cells_in_region * 100)

astro_region_data$broad_region <- factor(astro_region_data$broad_region, levels = c("Subcortical", "Cortical"))

astro_region_plot <- ggplot(astro_region_data, aes(x = prop, y = cluster_anno, fill = broad_region)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  scale_fill_manual(values = c(Cortical = "#9e1f63", Subcortical = "#6fc251")) +
  scale_y_discrete(limits = astro_cluster_order) +
  scale_x_continuous(position = "top", transform = "log1p") +
  theme_empty() +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 10, angle = 60, hjust = 0), 
        axis.ticks.x = element_line(linewidth = 0.2), axis.line.x = element_line(linewidth = 0.2),
        text = element_text(family = "Arial")) 

##################################################

# dotplot of marker genes

astro_genes <- unique(c("AQP4", "ALDH1L1", "ADGRV1", "GFAP", "S100B", 
                        "CD44", "SLC38A1", "CABLES1", "SLC1A2", 
                        "GRM3", "ARHGAP24", "WIF1", "RERG", "PDE3B", "VAV3",
                        "DCLK1", "CCDC85A", "ROBO2", "L3MBTL4", "ADAMTSL3", "NFASC", "SYNPO2",
                        "GRIA1", "LINC01411", "SEMA5A", "GUCY1A1", "KCNN2",
                        "KCND2", "LUZP2", "RALYL", "STK32A", "KIAA1217", "GRID2", "LGR6",
                        "SERPINA3", "CHI3L1", "PLSCR1", "OSMR", "SOCS3", "ZFP36", "NDRG1", "ATP11A",
                        "LAMA2", "STXBP5L", "GRIN2A", "SPAG17", "DNAH11",
                        "HSPH1", "HSPA6", "DNAJA4", "MMP16", "SERPINH1", "UBB", "CRYAB"))

astro_dotplot <- DotPlot(astro, features = astro_genes, assay = "sketch", dot.min = 0.05, dot.scale = 6) +
  aes(colour = avg.exp.scaled) + 
  scale_color_gradient(limits = c(-1, 2), low = "#d9f4fa", high = "#64288a", oob = squish) +
  scale_x_discrete(limits = astro_genes, position = "top") +
  scale_y_discrete(limits = astro_cluster_order) +
  theme(legend.position = "right", axis.title = element_blank(), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.2), axis.ticks = element_line(linewidth = 0.2), 
        axis.text.x = element_text(size = 10, angle = 60, hjust = 0), axis.text.y = element_blank(),
        text = element_text(family = "Arial")) +
  NoLegend()

########

# combine plots

p2 <- plot_grid(astro_ncells_plot, astro_region_plot, astro_dotplot, nrow = 1, ncol = 3, rel_widths = c(1.5, 1.5, 7),
                align = "h", axis = "tb") %>% print()

##################################################
##################################################
##################################################

####################
# OLIGODENDROCYTES #
####################

oligos <- readRDS("./analysis/oligodendrocytes/seurat_objects/oligos_clustered_ANNOTATED.rds")
oligos_meta <- readRDS("./analysis/oligodendrocytes/oligos_celllevel_metadata_ANNOTATED.rds")

oligos_cluster_order <- rev(c("Oligo_OPALIN", "Oligo_RBFOX1", "Oligo_CSMD1", "Oligo_CNTN1", "Oligo_Stress_HSPH1"))

##################################################

# number of cells, plotted by dataset

oligos_ncells_data <- oligos_meta %>%
  group_by(cluster_anno, dataset) %>%
  summarize(ncells = n()) %>%
  ungroup()

oligos_ncells_data$dataset <- ifelse(oligos_ncells_data$dataset == "Gerrits", paste0(oligos_ncells_data$dataset, "_2022"),
                                    ifelse(oligos_ncells_data$dataset == "AMP-PD", "NM_2024", paste0(oligos_ncells_data$dataset, "_2024")))

oligos_ncells_plot <- ggplot(oligos_ncells_data, aes(x = ncells, y = cluster_anno, fill = dataset)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = study_colors) +
  scale_y_discrete(limits = oligos_cluster_order) +
  scale_x_continuous(position = "top", labels = label_comma()) +
  theme_empty() +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 10, angle = 60, hjust = 0), 
        axis.ticks.x = element_line(linewidth = 0.2), 
        axis.line.x = element_line(linewidth = 0.2),
        axis.line.y = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2), 
        text = element_text(family = "Arial"))

##################################################

# proportions of cells cortical/subcortical

oligos_meta$broad_region <- ifelse(oligos_meta$region %in% c("DMV", "GPi", "TH"), "Subcortical", "Cortical")

oligos_region_data <- oligos_meta %>%
  group_by(cluster_anno, broad_region) %>%
  summarise(n = n()) %>% 
  group_by(broad_region) %>%
  mutate(num_cells_in_region = sum(n)) %>%
  mutate(prop = n / num_cells_in_region * 100)

oligos_region_data$broad_region <- factor(oligos_region_data$broad_region, levels = c("Subcortical", "Cortical"))

oligos_region_plot <- ggplot(oligos_region_data, aes(x = prop, y = cluster_anno, fill = broad_region)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  scale_fill_manual(values = c(Cortical = "#9e1f63", Subcortical = "#6fc251")) +
  scale_y_discrete(limits = oligos_cluster_order) +
  scale_x_continuous(position = "top", transform = "log1p") +
  theme_empty() +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 10, angle = 60, hjust = 0), 
        axis.ticks.x = element_line(linewidth = 0.2), axis.line.x = element_line(linewidth = 0.2),
        text = element_text(family = "Arial")) 

##################################################

# dotplot of marker genes

oligos_genes <- unique(c("OPALIN", "LAMA2", "PLXDC2", "DYSF", "LINC01608", "FCHSD2",
                         "RBFOX1", "RASGRF1", "AFF3", "FSTL5", "KANK4", "PCLO", "MBOAT1", "TNFRSF21",
                         "CSMD1", "SPARCL1", "MEG3", "SYT1", "KCNIP4", "DPP10", "NRXN1",
                         "CNTN1", "LUZP2", "SGCZ", "GRIK2", "FRY", "MDGA2", "ZFPM2", "CSMD3",
                         "HSPH1", "HSPA1B", "FOS", "HSPA1B", "BAG3"))

oligos_dotplot <- DotPlot(oligos, features = oligos_genes, assay = "sketch", dot.min = 0.05, dot.scale = 6) +
  aes(colour = avg.exp.scaled) +
  scale_color_gradient(limits = c(-1, 2), low = "#d9f4fa", high = "#64288a", oob = squish) +
  scale_x_discrete(limits = oligos_genes, position = "top") +
  scale_y_discrete(limits = oligos_cluster_order) +
  theme(legend.position = "right", axis.title = element_blank(), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.2), axis.ticks = element_line(linewidth = 0.2), 
        axis.text.x = element_text(size = 10, angle = 60, hjust = 0), axis.text.y = element_blank(),
        text = element_text(family = "Arial")) +
  NoLegend()

########

# combine plots

p3 <- plot_grid(oligos_ncells_plot, oligos_region_plot, oligos_dotplot, nrow = 1, ncol = 3, rel_widths = c(1.5, 1.5, 7),
                align = "h", axis = "tb") %>% print()

##################################################
##################################################
##################################################

# merge all plots

p4 <- plot_grid(p2, p1, p3, ncol = 1, rel_heights = c(6.5, 10, 5))

svglite(filename = "./analysis/misc/celltype_summary_plot_full.svg", width = 15.5, height = 9)
plot(p4)
dev.off()

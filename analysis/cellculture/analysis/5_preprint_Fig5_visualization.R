library(tidyverse)
library(fgsea)
library(eulerr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(ComplexHeatmap)
library(circlize)
library(purrr)
library(patchwork)
library(svglite)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")

##################################################
##################################################
##################################################

XR10585 <- readRDS("./analysis/cellculture/iMG_bulk_RNA/tmp_results/XR10585_dream_res_df_list.rds")
DA10634 <- readRDS("./analysis/cellculture/iMG_bulk_RNA/tmp_results/DA10634_dream_res_df_list.rds")


dream_res_list <- c(XR10585, DA10634)

##################################################
##################################################
##################################################

# filter all dfs for DAM genes + make heatmap

DAM_res <- list()
logFC_list <- list()
fdr_list <- list()

for (test in names(dream_res_list)){
  df <- dream_res_list[[test]]
  
  DAM_df <- df[rownames(df) %in% DAM_genes, ]
  
  DAM_res[[test]] <- DAM_df
  
  
  df1 <- DAM_df %>%
    rownames_to_column(var = "gene") %>%
    dplyr::select(gene, logFC) %>%
    dplyr::rename(!!paste0(test) := logFC)
  
  df2 <- DAM_df %>%
    rownames_to_column(var = "gene") %>%
    dplyr::select(gene, adj.P.Val) %>%
    dplyr::rename(!!paste0(test) := adj.P.Val)
  
  logFC_list[[test]] <- df1
  fdr_list[[test]] <- df2
}

#########################

DAM_res_df <- purrr::reduce(logFC_list, full_join, by = "gene")

DAM_res_df <- DAM_res_df %>%
  column_to_rownames(var = "gene") %>%
  dplyr::select(-`Vorinostat_1e-05`, -`Entinostat_1e-08`) %>%
  as.matrix()

DAM_res_df <- DAM_res_df[, grepl("UV-DAMPS|Entinostat|Vorinostat", colnames(DAM_res_df))]

#########################

DAM_fdr_df <- purrr::reduce(fdr_list, full_join, by = "gene")

DAM_fdr_df <- DAM_fdr_df %>%
  column_to_rownames(var = "gene") %>%
  dplyr::select(-`Vorinostat_1e-05`, -`Entinostat_1e-08`) %>%
  as.matrix()

DAM_fdr_df <- DAM_fdr_df[, grepl("UV-DAMPS|Entinostat|Vorinostat", colnames(DAM_fdr_df))]

##################################################
##################################################
##################################################

# run FGSEA for hnDAM/apoptosis sigs

# apoptosis: hsa04210

pathway_id <- "hsa04210"

# retrieve pathway data
pathway_data <- keggGet(pathway_id)

# extract Entrez Gene IDs
genes_raw <- pathway_data[[1]]$GENE
entrez_ids <- genes_raw[seq(1, length(genes_raw), 2)]

# convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

fgsea_res <- data.frame()

set.seed(12345)

for (test in names(dream_res_list)){
  test_df <- dream_res_list[[test]]
  
  list <- -log10(test_df$P.Value) * sign(test_df$logFC)
  names(list) <- rownames(test_df)
  
  fgseares <- fgsea(pathways = list("apoptosis" = gene_symbols,
                                    "hnDAM" = DAM_genes),
                    stats = list)
  
  fgsea_res <- rbind(fgsea_res, 
                     data.frame(test = test,
                                NES = fgseares$NES,
                                pval = fgseares$pval,
                                pathway = fgseares$pathway))
}

fgsea_res <- fgsea_res[fgsea_res$test %in% colnames(DAM_res_df), ]

fgsea_res$FDR <- p.adjust(fgsea_res$pval, method = "BH")

write.csv(fgsea_res, "./analysis/csvs_for_supp_tables/iMG_fgsea.csv", row.names = F)

##################################################
##################################################
##################################################

# make heatmap

##################################################

col_fun <- colorRamp2(breaks = c(-2, 0, 2), colors = c("#05409e", "white", "#e64a02"))
col_fun2 <- colorRamp2(breaks = c(-2, 0, 2), colors = c("#a6611a", "white", "#018571"))

##################################################

# make heatmap annotation for fgsea results

anno_df <- fgsea_res %>%
  dplyr::select(-pval, -FDR) %>%
  pivot_wider(names_from = pathway, values_from = NES)

anno_fdr_df <- fgsea_res %>%
  dplyr::select(-pval, -NES) %>%
  pivot_wider(names_from = pathway, values_from = FDR)


# function to overlay significance stars for GSEA results
annotation_with_stars <- function(values, fdr_values, col_fun) {
  anno_simple(
    values,
    col = col_fun,
    pch = ifelse(fdr_values <= 0.05, "﹡", ""),
    pt_gp = gpar(fontsize = 12, col = "black"),
    gp = gpar(col = "black"), 
    border = TRUE           
  )
}


ra <- rowAnnotation(
  hnDAM = annotation_with_stars(anno_df$hnDAM, anno_fdr_df$hnDAM, col_fun2),
  Apoptosis = annotation_with_stars(anno_df$apoptosis, anno_fdr_df$apoptosis, col_fun2),
  show_legend = FALSE)


# custom legend for NES
nes_legend <- Legend(col_fun = col_fun2,
                     title = "NES")

##################################################

fdr_fun <- function(j, i, x, y, w, h, fill) {
  if(t(DAM_fdr_df)[i, j] <= 0.05) {
    grid.text("﹡", x, y)
  }
}

##################################################

colnames(DAM_res_df) <- gsub("1e-07", "0.1μM", colnames(DAM_res_df))
colnames(DAM_res_df) <- gsub("1e-06", "1μM", colnames(DAM_res_df))

##################################################

ht = Heatmap(t(DAM_res_df), 
             cluster_rows = F, 
             col = col_fun, 
             right_annotation = ra,
             cell_fun = fdr_fun,
             name = "logFC")


png("./analysis/cellculture/iMG_bulk_RNA/plots_final/DAM_genes_heatmap.png", width = 20, height = 5, units = "in", res = 600)
draw(ht, annotation_legend_list = list(nes_legend), merge_legend = T)
dev.off()

##################################################
##################################################
##################################################

# barplot for number of sig DEGs

nsigs_df <- data.frame()

for (test in names(dream_res_list)){
  df <- dream_res_list[[test]]
  
  n_sigs <- as.numeric(nrow(df[df$adj.P.Val < 0.05 & abs(df$logFC) > 0.5, ]))
  
  nsigs_df <- rbind(nsigs_df, 
                    data.frame(test = test,
                               n_sigs = n_sigs))
}

nsigs_df$test <- gsub("1e-07", "0.1μM", nsigs_df$test)
nsigs_df$test <- gsub("1e-06", "1μM", nsigs_df$test)


nsigs_df <- nsigs_df[nsigs_df$test %in% colnames(DAM_res_df), ]

nsigs_df$test <- factor(nsigs_df$test, levels = rev(nsigs_df$test))


p1 <- ggplot(nsigs_df, aes(x = n_sigs, y = test, fill = test)) +
  geom_bar(stat = "identity")+
  theme_bw() +
  scale_fill_manual(values = rev(c("#edcccc", "#d6594b", "#821717",
                                   "#daedcc", "#2e8217",
                                   "#ccd9ed", "#172c82"))) + 
  labs(x = "# DEGs (FDR < 0.05, |logFC| > 0.5)",
       y = "Exposure") +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "none")

svglite("./analysis/cellculture/iMG_bulk_RNA/plots_final/nsigs_barplot.svg", width = 5.6, height = 3.8)
plot(p1)
dev.off()

##################################################
##################################################
##################################################

# making GSEA plots -- combining DAMPS and HDACs onto one plot

##################################################
##################################################
##################################################

# get plotting data for all conditions

curve_data_list <- list()
ticks_data_list <- list()

pathway_list <- list(hnDAM = DAM_genes,
                     apoptosis = gene_symbols)

for (test in names(dream_res_list)[c(1:3, 15:18)]){
  test_df <- dream_res_list[[test]]
  list <- -log10(test_df$P.Value) * sign(test_df$logFC)
  names(list) <- rownames(test_df)
  
  # loop thru pathways to extract plot data
  for (pathway in names(pathway_list)){
    test_pathway <- paste0(test, "_", pathway)
    
    path <- as.character(pathway_list[[pathway]])
    
    # get enrichment plot data
    plot_data <- fgsea::plotEnrichmentData(pathway = path, stats = list)
    
    # running enrichment score curve data
    curve <- plot_data$curve %>%
      as.data.frame() %>%
      dplyr::select(rank = 1, value = 2) %>%
      mutate(curve = test,
             pathway = pathway)
    curve_data_list[[test_pathway]] <- curve
    
    # gene ticks data
    ticks <- data.frame(rank = plot_data$ticks$rank,
                        curve = test,
                        pathway = pathway)
    ticks_data_list[[test_pathway]] <- ticks
  }
}

curve_data <- bind_rows(curve_data_list)
ticks_data <- bind_rows(ticks_data_list)

##################################################
##################################################
##################################################

# plotting curves and ticks

##################################################

# DAMPS -- hnDAM

DAMPS_hnDAM_curve <- curve_data %>%
  filter(grepl("UV-DAMPS", curve) & pathway == "hnDAM")
DAMPS_hnDAM_ticks <- ticks_data %>%
  filter(grepl("UV-DAMPS", curve) & pathway == "hnDAM")

condition_levels <- c("UV-DAMPS_168000", "UV-DAMPS_336000", "UV-DAMPS_672000")
condition_colors <- c("UV-DAMPS_168000" = "#edcccc", 
                      "UV-DAMPS_336000" = "#d6594b", 
                      "UV-DAMPS_672000" = "#821717")
names(condition_colors) <- condition_levels


DAMPS_hnDAM_curve$curve <- factor(DAMPS_hnDAM_curve$curve, levels = condition_levels)
DAMPS_hnDAM_ticks$curve <- factor(DAMPS_hnDAM_ticks$curve, levels = condition_levels)


DAMPS_hnDAM_curve_plot <- ggplot(DAMPS_hnDAM_curve, aes(x = rank, y = value, color = curve)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(y = "Running enrichment score", 
       title = "UV-DAMPS - hnDAM") +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(limits = c(-0.45, 0.45)) + 
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "none")


DAMPS_hnDAM_ticks_plot <- ggplot(DAMPS_hnDAM_ticks, aes(x = rank, y = curve, color = curve)) +
  geom_segment(aes(xend = rank, y = as.numeric(curve) - 0.4, yend = as.numeric(curve) + 0.4), linewidth = 0.5) +
  scale_color_manual(values = condition_colors) +
  labs(x = "Rank",
       color = "Exposure") +
  theme_bw() + 
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 1),
    legend.position = "bottom",
    legend.direction = "horizontal")


DAMPS_hnDAM_merged_plot <- (DAMPS_hnDAM_curve_plot / DAMPS_hnDAM_ticks_plot) +
  plot_layout(heights = c(4, 1), guides = 'collect')

DAMPS_hnDAM_merged_plot

##################################################

# DAMPS -- apoptosis

DAMPS_apoptosis_curve <- curve_data %>%
  filter(grepl("UV-DAMPS", curve) & pathway == "apoptosis")
DAMPS_apoptosis_ticks <- ticks_data %>%
  filter(grepl("UV-DAMPS", curve) & pathway == "apoptosis")

condition_levels <- c("UV-DAMPS_168000", "UV-DAMPS_336000", "UV-DAMPS_672000")
condition_colors <- c("UV-DAMPS_168000" = "#edcccc", 
                      "UV-DAMPS_336000" = "#d6594b", 
                      "UV-DAMPS_672000" = "#821717")
names(condition_colors) <- condition_levels


DAMPS_apoptosis_curve$curve <- factor(DAMPS_apoptosis_curve$curve, levels = condition_levels)
DAMPS_apoptosis_ticks$curve <- factor(DAMPS_apoptosis_ticks$curve, levels = condition_levels)


DAMPS_apoptosis_curve_plot <- ggplot(DAMPS_apoptosis_curve, aes(x = rank, y = value, color = curve)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(y = "Running enrichment score", 
       title = "UV-DAMPS - Apoptosis (KEGG hsa04210)") +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(limits = c(-0.45, 0.45)) + 
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "none")


DAMPS_apoptosis_ticks_plot <- ggplot(DAMPS_apoptosis_ticks, aes(x = rank, y = curve, color = curve)) +
  geom_segment(aes(xend = rank, y = as.numeric(curve) - 0.4, yend = as.numeric(curve) + 0.4), linewidth = 0.5) +
  scale_color_manual(values = condition_colors) +
  labs(x = "Rank",
       color = "Exposure") +
  theme_bw() + 
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 1),
    legend.position = "bottom",
    legend.direction = "horizontal")


DAMPS_apoptosis_merged_plot <- DAMPS_apoptosis_curve_plot / DAMPS_apoptosis_ticks_plot +
  plot_layout(heights = c(4, 1), guides = 'collect')

DAMPS_apoptosis_merged_plot

##################################################

# HDAC -- hnDAM

HDAC_hnDAM_curve <- curve_data %>%
  filter(grepl("Entinostat|Vorinostat", curve) & pathway == "hnDAM")
HDAC_hnDAM_ticks <- ticks_data %>%
  filter(grepl("Entinostat|Vorinostat", curve) & pathway == "hnDAM")


HDAC_hnDAM_curve$curve <- gsub("1e-07", "0.1μM", HDAC_hnDAM_curve$curve)
HDAC_hnDAM_curve$curve <- gsub("1e-06", "1μM", HDAC_hnDAM_curve$curve)
HDAC_hnDAM_ticks$curve <- gsub("1e-07", "0.1μM", HDAC_hnDAM_ticks$curve)
HDAC_hnDAM_ticks$curve <- gsub("1e-06", "1μM", HDAC_hnDAM_ticks$curve)

condition_levels <- c("Entinostat_0.1μM", "Entinostat_1μM", "Vorinostat_0.1μM", "Vorinostat_1μM")
condition_colors <- c("Entinostat_0.1μM" = "#daedcc", 
                      "Entinostat_1μM" = "#2e8217", 
                      "Vorinostat_0.1μM" = "#ccd9ed", 
                      "Vorinostat_1μM" = "#172c82")
names(condition_colors) <- condition_levels


HDAC_hnDAM_curve$curve <- factor(HDAC_hnDAM_curve$curve, levels = condition_levels)
HDAC_hnDAM_ticks$curve <- factor(HDAC_hnDAM_ticks$curve, levels = condition_levels)


HDAC_hnDAM_curve_plot <- ggplot(HDAC_hnDAM_curve, aes(x = rank, y = value, color = curve)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(y = "Running enrichment score", 
       title = "HDAC inhibitors - hnDAM") +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(limits = c(-0.45, 0.45)) + 
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "none")


HDAC_hnDAM_ticks_plot <- ggplot(HDAC_hnDAM_ticks, aes(x = rank, y = curve, color = curve)) +
  geom_segment(aes(xend = rank, y = as.numeric(curve) - 0.4, yend = as.numeric(curve) + 0.4), linewidth = 0.5) +
  scale_color_manual(values = condition_colors) +
  labs(x = "Rank",
       color = "Exposure") +
  theme_bw() + 
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 1),
    legend.position = "bottom",
    legend.direction = "horizontal")


HDAC_hnDAM_merged_plot <- HDAC_hnDAM_curve_plot / HDAC_hnDAM_ticks_plot +
  plot_layout(heights = c(4, 1), guides = 'collect')

HDAC_hnDAM_merged_plot

##################################################

# HDAC -- apoptosis

HDAC_apoptosis_curve <- curve_data %>%
  filter(grepl("Entinostat|Vorinostat", curve) & pathway == "apoptosis")
HDAC_apoptosis_ticks <- ticks_data %>%
  filter(grepl("Entinostat|Vorinostat", curve) & pathway == "apoptosis")


HDAC_apoptosis_curve$curve <- gsub("1e-07", "0.1μM", HDAC_apoptosis_curve$curve)
HDAC_apoptosis_curve$curve <- gsub("1e-06", "1μM", HDAC_apoptosis_curve$curve)
HDAC_apoptosis_ticks$curve <- gsub("1e-07", "0.1μM", HDAC_apoptosis_ticks$curve)
HDAC_apoptosis_ticks$curve <- gsub("1e-06", "1μM", HDAC_apoptosis_ticks$curve)

condition_levels <- c("Entinostat_0.1μM", "Entinostat_1μM", "Vorinostat_0.1μM", "Vorinostat_1μM")
condition_colors <- c("Entinostat_0.1μM" = "#daedcc", 
                      "Entinostat_1μM" = "#2e8217", 
                      "Vorinostat_0.1μM" = "#ccd9ed", 
                      "Vorinostat_1μM" = "#172c82")
names(condition_colors) <- condition_levels


HDAC_apoptosis_curve$curve <- factor(HDAC_apoptosis_curve$curve, levels = condition_levels)
HDAC_apoptosis_ticks$curve <- factor(HDAC_apoptosis_ticks$curve, levels = condition_levels)


HDAC_apoptosis_curve_plot <- ggplot(HDAC_apoptosis_curve, aes(x = rank, y = value, color = curve)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(y = "Running enrichment score", 
       title = "HDAC inhibitors - Apoptosis (KEGG hsa04210)") +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(limits = c(-0.45, 0.45)) + 
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "none")


HDAC_apoptosis_ticks_plot <- ggplot(HDAC_apoptosis_ticks, aes(x = rank, y = curve, color = curve)) +
  geom_segment(aes(xend = rank, y = as.numeric(curve) - 0.4, yend = as.numeric(curve) + 0.4), linewidth = 0.5) +
  scale_color_manual(values = condition_colors) +
  labs(x = "Rank",
       color = "Exposure") +
  theme_bw() + 
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 1),
    legend.position = "bottom",
    legend.direction = "horizontal")


HDAC_apoptosis_merged_plot <- HDAC_apoptosis_curve_plot / HDAC_apoptosis_ticks_plot +
  plot_layout(heights = c(4, 1))

HDAC_apoptosis_merged_plot

##################################################

final_plot <- ((DAMPS_hnDAM_curve_plot / DAMPS_hnDAM_ticks_plot + 
                  plot_layout(heights = c(4, 1))) | 
                 (DAMPS_apoptosis_curve_plot / DAMPS_apoptosis_ticks_plot + 
                    plot_layout(heights = c(4, 1))) | 
                 (HDAC_hnDAM_curve_plot / HDAC_hnDAM_ticks_plot + 
                    plot_layout(heights = c(4, 1))) | 
                 (HDAC_apoptosis_curve_plot / HDAC_apoptosis_ticks_plot + 
                    plot_layout(heights = c(4, 1)))) + 
  plot_layout(ncol = 4)

final_plot

svglite("./analysis/cellculture/iMG_bulk_RNA/plots_final/iMG_gsea_combined_plot.svg", width = 20, height = 5)
plot(final_plot)
dev.off()

library(tidyverse)
library(fgsea)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# run on NON-ANNOTED CLUSTERS (to be used for annotation purposes)

##################################################

cortical <- readRDS("./analysis/astrocytes/seurat_objects/cortical_astrocytes_sketched_filtered_clustered_projected_0.15.rds")

cortical[["RNA"]] <- JoinLayers(cortical[["RNA"]])

cortical_markers_all <- FindAllMarkers(cortical, group.by = "celltype.full", assay = "RNA", only.pos = F, logfc.threshold = 0, min.pct = 0.2)
saveRDS(cortical_markers_all, file = "./analysis/astrocytes/cluster_characterization/FGSEA/cortical_markers_ALL_for_gsea.rds")

##################################################

astro_markers <- readRDS("./analysis/astrocytes/cluster_characterization/FGSEA/cortical_markers_ALL_for_gsea.rds")

##################################################

astro_signatures <- read.csv("./analysis/astrocytes/cluster_characterization/FGSEA/astro_signatures_to_test.csv", header = T)
astro_signatures_list <- lapply(astro_signatures, function(x) x[x != ""])
names(astro_signatures_list) <- names(astro_signatures)

##################################################

astro_clusters_list <- unique(astro_markers$cluster)

fgsea_res_list <- list()

for (i in astro_clusters_list){
  markers <- astro_markers[astro_markers$cluster == i, ]
  
  rankings <- sign(markers$avg_log2FC)*(-log10(markers$p_val_adj)) 
  names(rankings) <- markers$gene 
  max_ranking <- max(rankings[is.finite(rankings)])
  min_ranking <- min(rankings[is.finite(rankings)])
  rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
  rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
  rankings <- sort(rankings, decreasing = TRUE)
  
  fgseaRes <- fgsea(pathways = astro_signatures_list, 
                    stats = rankings, 
                    minSize = 0, 
                    maxSize = 500,
                    scoreType = "std")
  
  fgsea_res_list[[i]] <- fgseaRes
}

fgsea_pval_df <- Reduce(function(x, y) merge(x, y, by = "signatures", all = TRUE),
                        lapply(names(fgsea_res_list), function(name) {
                          df <- fgsea_res_list[[name]]
                          data.frame(signatures = df$pathway, value = df[[3]], check.names = FALSE)
                        }))

colnames(fgsea_pval_df)[-1] <- names(fgsea_res_list)
fgsea_pval_df <- fgsea_pval_df %>%
  column_to_rownames(var = "signatures")



fgsea_NES_df <- Reduce(function(x, y) merge(x, y, by = "signatures", all = TRUE),
                       lapply(names(fgsea_res_list), function(name) {
                         df <- fgsea_res_list[[name]]
                         data.frame(signatures = df$pathway, value = df[[6]], check.names = FALSE)
                       }))

colnames(fgsea_NES_df)[-1] <- names(fgsea_res_list)
fgsea_NES_df <- fgsea_NES_df %>%
  column_to_rownames(var = "signatures")


fgsea_pval_df <- fgsea_pval_df %>%
  rownames_to_column(var = "signature") %>%
  pivot_longer(cols = -signature, names_to = "cluster", values_to = "pval")

fgsea_NES_df <- fgsea_NES_df %>%
  rownames_to_column(var = "signature") %>%
  pivot_longer(cols = -signature, names_to = "cluster", values_to = "NES")

fgsea_df <- left_join(fgsea_pval_df, fgsea_NES_df, by = c("cluster", "signature"))

fgsea_df$padj <- p.adjust(fgsea_df$pval, method = "BH")


fgsea_df$log10padj <- -log10(fgsea_df$padj)
fgsea_df$significant <- ifelse(fgsea_df$padj <= 0.05, T, F)

saveRDS(fgsea_df, "./analysis/astrocytes/cluster_characterization/FGSEA/cortical_astro_FGSEA_NOT_annotated.rds")

##################################################

# visualize for non-annotated clusters

fgsea_df <- readRDS("./analysis/astrocytes/cluster_characterization/FGSEA/cortical_astro_FGSEA_NOT_annotated.rds")

ggplot(fgsea_df, aes(x = cluster, y = signature)) +
  geom_point(aes(size = log10padj, fill = NES, color = significant), shape = 21, stroke = 0.5) +
  scale_size(range = c(1, 10)) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") + 
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Cluster", y = "Signature", size = "-log10(padj)", color = "Significant", fill = "Normalized Enrichment") +
  coord_flip()


library(tidyverse)
library(fgsea)
library(Seurat)
library(ggplot2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# run on NON-ANNOTED CLUSTERS (to be used for annotation purposes)

##################################################

micro <- readRDS("./combined_data/final_objects/micro_sketched_clustered_projected_NO_ctsfilter_res_0.25.rds")
micro[["RNA"]] <- JoinLayers(micro[["RNA"]])

micro_markers_all <- FindAllMarkers(micro, group.by = "celltype.full", assay = "RNA", only.pos = F,
                                    logfc.threshold = 0, min.pct = 0.2)
saveRDS(micro_markers_all, file = "./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea.rds")
write.csv(micro_markers_all, file = "./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea.csv")

##################################################

micro_markers <- readRDS("./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea.rds")

##################################################

micro_signatures <- read.csv("./analysis/microglia/cluster_characterization/FGSEA/micro_signatures_to_test.csv", header = T)
micro_signatures <- micro_signatures[, -grep("Drager", colnames(micro_signatures))]
micro_signatures_list <- lapply(micro_signatures, function(x) x[x != ""])
names(micro_signatures_list) <- names(micro_signatures)

##################################################

micro_clusters_list <- c("0", "1", "2", "3", "4", "5", "6", "7", 
                         "8", "9", "10", "11", "12", "13", "14")

fgsea_res_list <- list()

for (i in micro_clusters_list){
  markers <- micro_markers[micro_markers$cluster == i, ]
  
  rankings <- sign(markers$avg_log2FC)*(-log10(markers$p_val_adj)) 
  names(rankings) <- markers$gene 
  max_ranking <- max(rankings[is.finite(rankings)])
  min_ranking <- min(rankings[is.finite(rankings)])
  rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
  rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
  rankings <- sort(rankings, decreasing = TRUE)
  
  fgseaRes <- fgsea(pathways = micro_signatures_list, 
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

##################################################

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
fgsea_df$cluster <- factor(fgsea_df$cluster, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                                        "10", "11", "12", "13", "14"))

saveRDS(fgsea_df, "./analysis/microglia/cluster_characterization/FGSEA/micro_FGSEA_signatures_output.rds")

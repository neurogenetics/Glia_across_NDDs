library(tidyverse)
library(fgsea)
library(Seurat)
library(ggplot2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

micro_markers <- readRDS("./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea_ANNOTATED.rds")

##################################################

# testing ALL signatures to see which ones to plot

micro_signatures <- read.csv("./analysis/microglia/cluster_characterization/FGSEA/micro_ALL_FGSEA_signatures_test.csv", header = T)
micro_signatures_list <- lapply(micro_signatures, function(x) x[x != ""])
names(micro_signatures_list) <- names(micro_signatures)

micro_signatures_list$Gerrits_2021_c7_AD1 <- gsub("\\.[^.]*$", "", micro_signatures_list$Gerrits_2021_c7_AD1)
micro_signatures_list$Gerrits_2021_c9_AD1 <- gsub("\\.[^.]*$", "", micro_signatures_list$Gerrits_2021_c9_AD1)
micro_signatures_list$Gerrits_2021_c10_AD1 <- gsub("\\.[^.]*$", "", micro_signatures_list$Gerrits_2021_c10_AD1)

# subset each signature to only top 100 values
# micro_signatures_list <- lapply(micro_signatures_list, function(x) x[1:min(100, length(x))])

##################################################

micro_clusters_list <- unique(as.character(micro_markers$cluster))

fgsea_res_list <- list()

# i = "Micro_DAM_GPNMB"

for (i in micro_clusters_list){
  markers <- micro_markers[micro_markers$cluster == i, ]
  
  markers$p_val_adj[markers$p_val_adj == 0] <- .Machine$double.xmin
  
  rankings <- markers$avg_log2FC*(-log10(markers$p_val_adj)) 
  names(rankings) <- markers$gene 
  # max_ranking <- max(rankings[is.finite(rankings)])
  # min_ranking <- min(rankings[is.finite(rankings)])
  # rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
  # rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
  rankings <- sort(rankings, decreasing = TRUE)
  
  fgseaRes <- fgsea(pathways = micro_signatures_list, 
                    stats = rankings, 
                    minSize = 0, 
                    maxSize = 1000,
                    scoreType = "std",
                    eps = 0,
                    nPermSimple = 10000)
  
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

##################################################
##################################################
##################################################

# run w/ final signatures to plot 

micro_markers <- readRDS("./analysis/microglia/cluster_characterization/FGSEA/micro_markers_ALL_for_gsea_ANNOTATED.rds")
micro_signatures <- read.csv("./analysis/microglia/cluster_characterization/FGSEA/micro_FINAL_signatures_to_test.csv", header = T)

##################################################

micro_signatures_list <- lapply(micro_signatures, function(x) x[x != ""])
names(micro_signatures_list) <- names(micro_signatures)

micro_signatures_list$Gerrits_2021_c7_AD1 <- gsub("\\.[^.]*$", "", micro_signatures_list$Gerrits_2021_c7_AD1)
micro_signatures_list$Gerrits_2021_c9_AD1 <- gsub("\\.[^.]*$", "", micro_signatures_list$Gerrits_2021_c9_AD1)
micro_signatures_list$Gerrits_2021_c10_AD1 <- gsub("\\.[^.]*$", "", micro_signatures_list$Gerrits_2021_c10_AD1)

##################################################

micro_clusters_list <- unique(as.character(micro_markers$cluster))
micro_clusters_list <- micro_clusters_list[grepl("Micro", micro_clusters_list)]

fgsea_res_list <- list()

# i = "Micro_DAM_GPNMB"

for (i in micro_clusters_list){
  markers <- micro_markers[micro_markers$cluster == i, ]
  
  markers$p_val_adj[markers$p_val_adj == 0] <- .Machine$double.xmin
  
  rankings <- markers$avg_log2FC*(-log10(markers$p_val_adj)) 
  names(rankings) <- markers$gene 
  rankings <- sort(rankings, decreasing = TRUE)
  
  fgseaRes <- fgsea(pathways = micro_signatures_list, 
                    stats = rankings, 
                    minSize = 0, 
                    maxSize = 1000,
                    scoreType = "std",
                    eps = 0,
                    nPermSimple = 10000)
  
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

#########################

saveRDS(fgsea_df, "./analysis/microglia/cluster_characterization/FGSEA/micro_FGSEA_signatures_output_ANNOTATED.rds")

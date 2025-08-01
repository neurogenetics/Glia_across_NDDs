library(tidyverse)
library(fgsea)
library(ggplot2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# CORTICAL

##################################################

cortical_markers <- readRDS("./analysis/astrocytes/cluster_characterization/cortical_astros_markers_ALL_ANNOTATED.rds")

##################################################

astro_signatures <- read.csv("./analysis/astrocytes/cluster_characterization/FGSEA/astro_signatures_to_test.csv", header = T)
astro_signatures_list <- lapply(astro_signatures, function(x) x[x != ""])
names(astro_signatures_list) <- names(astro_signatures)

##################################################

astro_clusters_list <- unique(as.character(cortical_markers$cluster))

fgsea_res_list <- list()

for (i in astro_clusters_list){
  markers <- cortical_markers[cortical_markers$cluster == i, ]
  
  markers$p_val_adj[markers$p_val_adj == 0] <- .Machine$double.xmin
  
  rankings <- markers$avg_log2FC*(-log10(markers$p_val_adj)) 
  names(rankings) <- markers$gene 
  rankings <- sort(rankings, decreasing = TRUE)
  
  fgseaRes <- fgsea(pathways = astro_signatures_list, 
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

saveRDS(fgsea_df, "./analysis/astrocytes/cluster_characterization/FGSEA/astro_FGSEA_signatures_output_ANNOTATED.rds")

##################################################
##################################################
##################################################

# ALL ASTROCYTES

##################################################

cortical_markers <- readRDS("./analysis/astrocytes/cluster_characterization/all_astros_markers_ALL_ANNOTATED.rds")

##################################################

astro_signatures <- read.csv("./analysis/astrocytes/cluster_characterization/FGSEA/astro_signatures_to_test.csv", header = T)
astro_signatures_list <- lapply(astro_signatures, function(x) x[x != ""])
names(astro_signatures_list) <- names(astro_signatures)

##################################################

astro_clusters_list <- unique(as.character(cortical_markers$cluster))

fgsea_res_list <- list()

for (i in astro_clusters_list){
  markers <- cortical_markers[cortical_markers$cluster == i, ]
  
  markers$p_val_adj[markers$p_val_adj == 0] <- .Machine$double.xmin
  
  rankings <- markers$avg_log2FC*(-log10(markers$p_val_adj)) 
  names(rankings) <- markers$gene 
  rankings <- sort(rankings, decreasing = TRUE)
  
  fgseaRes <- fgsea(pathways = astro_signatures_list, 
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

saveRDS(fgsea_df, "./analysis/astrocytes/cluster_characterization/FGSEA/ALL_astros_FGSEA_signatures_output_ANNOTATED.rds")

library(corrplot)
library(ComplexHeatmap)
library(tidyverse)
library(mashr)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

setwd("/data/ADRD/glia_across_NDDs")

########################################

astro_mashr <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_diffexp_mashr_object.rds")
strong_subset_genes <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_strong_subset_genes.rds")

########################################

# get genes that fit different covariance matrices the best

postweights <- as.data.frame(astro_mashr$posterior_weights)
rownames(postweights) <- rownames(astro_mashr$result$PosteriorMean)
postweights <- postweights %>%
  rownames_to_column(var = "gene")

best_fit <- data.frame(gene = postweights$gene,
                       best_fit = colnames(postweights[, -1])[max.col(postweights[, -1], ties.method = "first")])


all_betas <- as.data.frame(astro_mashr$result$PosteriorMean)

all_lfsrs <- as.data.frame(astro_mashr$result$lfsr)

# all lfsrs for genes in strong subset
all_lfsrs_strong <- all_lfsrs[rownames(all_lfsrs) %in% strong_subset_genes, ]

########################################

# tPCA genes

tPCA_bestfit <- best_fit$gene[grepl("tPCA", best_fit$best_fit)]
tPCA_betas <- all_betas[rownames(all_betas) %in% tPCA_bestfit, ]

tPCA_betas <- tPCA_betas[rownames(tPCA_betas) %in% rownames(all_lfsrs_strong), ]

tPCA_lfsrs <- all_lfsrs_strong[rownames(all_lfsrs_strong) %in% rownames(tPCA_betas), ]

####################

saveRDS(tPCA_betas, "./analysis/astrocytes/differential_expression/astrocytes_mashr_tPCA_betas.rds")
saveRDS(tPCA_lfsrs, "./analysis/astrocytes/differential_expression/astrocytes_mashr_tPCA_lfsrs.rds")

########################################

# ED_PCA_1 genes

PCA1_bestfit <- best_fit$gene[grepl("ED_PCA_1", best_fit$best_fit)]
PCA1_betas <- all_betas[rownames(all_betas) %in% PCA1_bestfit, ]

PCA1_betas <- PCA1_betas[rownames(PCA1_betas) %in% rownames(all_lfsrs_strong), ]

PCA1_lfsrs <- all_lfsrs_strong[rownames(all_lfsrs_strong) %in% rownames(PCA1_betas), ]

####################

saveRDS(PCA1_betas, "./analysis/astrocytes/differential_expression/astrocytes_mashr_PCA1_betas.rds")
saveRDS(PCA1_lfsrs, "./analysis/astrocytes/differential_expression/astrocytes_mashr_PCA1_lfsrs.rds")

########################################

# ED_PCA_2 genes

PCA2_bestfit <- best_fit$gene[grepl("ED_PCA_2", best_fit$best_fit)]
PCA2_betas <- all_betas[rownames(all_betas) %in% PCA2_bestfit, ]

PCA2_betas <- PCA2_betas[rownames(PCA2_betas) %in% rownames(all_lfsrs_strong), ]

PCA2_lfsrs <- all_lfsrs_strong[rownames(all_lfsrs_strong) %in% rownames(PCA2_betas), ]

####################

saveRDS(PCA2_betas, "./analysis/astrocytes/differential_expression/astrocytes_mashr_PCA2_betas.rds")
saveRDS(PCA2_lfsrs, "./analysis/astrocytes/differential_expression/astrocytes_mashr_PCA2_lfsrs.rds")

########################################
########################################
########################################

# make simple heatmap and get k-means clusters of genes

tPCA_betas <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_mashr_tPCA_betas.rds")

# determine optimal K for k-means clustering

wss <- sapply(1:10, function(k){
  kmeans(tPCA_betas, centers = k, nstart = 25)$tot.withinss
})

plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster Sum of Squares")


# make and initialize heatmap

quantile(as.matrix(tPCA_betas), probs=0.95) # 1.264841
quantile(as.matrix(tPCA_betas), probs=0.05) # -1.319144

col_fun <- colorRamp2(c(-1, 0, 1), c("#05409e", "white", "#e64a02"))

set.seed(12345)
hmp1 <- Heatmap(as.matrix(tPCA_betas),
                show_row_names = F,
                cluster_rows = F,
                row_km = 6, 
                name = "Posterior mean",
                use_raster = F,
                col = col_fun)
hmp1 <- draw(hmp1)


# get genes from different k-clusters

for (i in 1:length(row_order(hmp1))){
  if (i == 1) {
    clu <- t(t(row.names(tPCA_betas[row_order(hmp1)[[i]],])))
    gene_clusters <- cbind(clu, paste("cluster", i, sep=""))
    colnames(gene_clusters) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(tPCA_betas[row_order(hmp1)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    gene_clusters <- rbind(gene_clusters, clu)
  }
}

gene_clusters <- as.data.frame(gene_clusters)

####################

saveRDS(gene_clusters, "./analysis/astrocytes/differential_expression/astrocytes_tPCA_kmeans_clusters.rds")

########################################
########################################
########################################

# run GO on k-means clusters

gene_clusters <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_tPCA_kmeans_clusters.rds")
astro_mashr <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_diffexp_mashr_object.rds")

####################

gobp_c2 <- enrichGO(gene = gene_clusters$GeneID[gene_clusters$Cluster == "cluster2"], 
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    keyType = "SYMBOL",
                    pAdjustMethod ="BH",
                    readable = TRUE,
                    universe = rownames(astro_mashr$result$lfsr))
gobp_c2 <- gobp_c2@result

gobp_c4 <- enrichGO(gene = gene_clusters$GeneID[gene_clusters$Cluster == "cluster4"], 
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    keyType = "SYMBOL",
                    pAdjustMethod ="BH",
                    readable = TRUE,
                    universe = rownames(astro_mashr$result$lfsr))
gobp_c4 <- gobp_c4@result

gobp_c5 <- enrichGO(gene = gene_clusters$GeneID[gene_clusters$Cluster == "cluster5"], 
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    keyType = "SYMBOL",
                    pAdjustMethod ="BH",
                    readable = TRUE,
                    universe = rownames(astro_mashr$result$lfsr))
gobp_c5 <- gobp_c5@result

####################

saveRDS(gobp_c2, "./analysis/astrocytes/differential_expression/astrocytes_gobp_c2.rds")
saveRDS(gobp_c4, "./analysis/astrocytes/differential_expression/astrocytes_gobp_c4.rds")
saveRDS(gobp_c5, "./analysis/astrocytes/differential_expression/astrocytes_gobp_c5.rds")

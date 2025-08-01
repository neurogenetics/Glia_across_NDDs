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
########################################
########################################

oligos_mashr <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_mashr_object.rds")
strong_subset_genes <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_strong_subset_genes.rds")

########################################
########################################
########################################

# get genes that fit different covariance matrices the best

postweights <- as.data.frame(oligos_mashr$posterior_weights)
rownames(postweights) <- rownames(oligos_mashr$result$PosteriorMean)
postweights <- postweights %>%
  rownames_to_column(var = "gene")

best_fit <- data.frame(gene = postweights$gene,
                       best_fit = colnames(postweights[, -1])[max.col(postweights[, -1], ties.method = "first")])


all_betas <- as.data.frame(oligos_mashr$result$PosteriorMean)

all_lfsrs <- as.data.frame(oligos_mashr$result$lfsr)

# all lfsrs for genes in strong subset
all_lfsrs_strong <- all_lfsrs[rownames(all_lfsrs) %in% strong_subset_genes, ]

########################################
########################################

# tPCA genes

tPCA_bestfit <- best_fit$gene[grepl("tPCA", best_fit$best_fit)]
tPCA_betas <- all_betas[rownames(all_betas) %in% tPCA_bestfit, ]

tPCA_betas <- tPCA_betas[rownames(tPCA_betas) %in% rownames(all_lfsrs_strong), ]

tPCA_lfsrs <- all_lfsrs_strong[rownames(all_lfsrs_strong) %in% rownames(tPCA_betas), ]

####################

saveRDS(tPCA_betas, "./analysis/oligodendrocytes/differential_expression/oligodendrocytes_mashr_tPCA_betas.rds")
saveRDS(tPCA_lfsrs, "./analysis/oligodendrocytes/differential_expression/oligodendrocytes_mashr_tPCA_lfsrs.rds")

########################################
########################################
########################################

# make simple heatmap and get k-means clusters of genes

tPCA_betas <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_mashr_tPCA_betas.rds")

# determine optimal K for k-means clustering

wss <- sapply(1:10, function(k){
  kmeans(tPCA_betas, centers = k, nstart = 25)$tot.withinss
})

plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster Sum of Squares")


# make and initialize heatmap

col_fun <- colorRamp2(c(-1, 0, 1), c("#05409e", "white", "#e64a02"))

set.seed(12345)
hmp1 <- Heatmap(as.matrix(tPCA_betas),
                show_row_names = F,
                cluster_rows = F,
                row_km = 5, 
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

saveRDS(gene_clusters, "./analysis/oligodendrocytes/differential_expression/oligodendrocytes_tPCA_kmeans_clusters.rds")

########################################
########################################
########################################

# run GO on k-means clusters

gene_clusters <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_tPCA_kmeans_clusters.rds")
oligo_mashr <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_mashr_object.rds")

####################

gobp_c2 <- enrichGO(gene = gene_clusters$GeneID[gene_clusters$Cluster == "cluster2"], 
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    keyType = "SYMBOL",
                    pAdjustMethod ="BH",
                    readable = TRUE,
                    universe = rownames(oligo_mashr$result$lfsr))
gobp_c2 <- gobp_c2@result



gobp_c3 <- enrichGO(gene = gene_clusters$GeneID[gene_clusters$Cluster == "cluster3"], 
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    keyType = "SYMBOL",
                    pAdjustMethod ="BH",
                    readable = TRUE,
                    universe = rownames(oligo_mashr$result$lfsr))
gobp_c3 <- gobp_c3@result


####################

saveRDS(gobp_c2, "./analysis/oligodendrocytes/differential_expression/oligodendrocytes_gobp_c2.rds")
saveRDS(gobp_c3, "./analysis/oligodendrocytes/differential_expression/oligodendrocytes_gobp_c3.rds")

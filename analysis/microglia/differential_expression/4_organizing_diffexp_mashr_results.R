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

micro_mashr <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_mashr_object.rds")
strong_subset_genes <- readRDS("./analysis/microglia/differential_expression/microglia_strong_subset_genes.rds")

########################################

# get genes that fit different covariance matrices the best

postweights <- as.data.frame(micro_mashr$posterior_weights)
rownames(postweights) <- rownames(micro_mashr$result$PosteriorMean)
postweights <- postweights %>%
  rownames_to_column(var = "gene")

best_fit <- data.frame(gene = postweights$gene,
                       best_fit = colnames(postweights[, -1])[max.col(postweights[, -1], ties.method = "first")])


all_betas <- as.data.frame(micro_mashr$result$PosteriorMean)

all_lfsrs <- as.data.frame(micro_mashr$result$lfsr)

# all lfsrs for genes in strong subset
all_lfsrs_strong <- all_lfsrs[rownames(all_lfsrs) %in% strong_subset_genes, ]

########################################

# tPCA genes

tPCA_bestfit <- best_fit$gene[grepl("tPCA", best_fit$best_fit)]
tPCA_betas <- all_betas[rownames(all_betas) %in% tPCA_bestfit, ]

tPCA_betas <- tPCA_betas[rownames(tPCA_betas) %in% rownames(all_lfsrs_strong), ]

tPCA_lfsrs <- all_lfsrs_strong[rownames(all_lfsrs_strong) %in% rownames(tPCA_betas), ]

####################

saveRDS(tPCA_betas, "./analysis/microglia/differential_expression/microglia_mashr_tPCA_betas.rds")
saveRDS(tPCA_lfsrs, "./analysis/microglia/differential_expression/microglia_mashr_tPCA_lfsrs.rds")

########################################
########################################
########################################

# make simple heatmap and get k-means clusters of genes

tPCA_betas <- readRDS("./analysis/microglia/differential_expression/microglia_mashr_tPCA_betas.rds")

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

saveRDS(gene_clusters, "./analysis/microglia/differential_expression/microglia_tPCA_kmeans_clusters.rds")

########################################
########################################
########################################

# run GO on k-means clusters

gene_clusters <- readRDS("./analysis/microglia/differential_expression/microglia_tPCA_kmeans_clusters.rds")
micro_mashr <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_mashr_object.rds")

####################

# run GO on k-means clusters

gobp_c1_c2 <- enrichGO(gene = gene_clusters$GeneID[gene_clusters$Cluster %in% c("cluster1", "cluster2")], 
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       keyType = "SYMBOL",
                       pAdjustMethod ="BH",
                       readable = TRUE,
                       universe = rownames(micro_mashr$result$lfsr))
gobp_c1_c2 <- gobp_c1_c2@result



gobp_c4_c6 <- enrichGO(gene = gene_clusters$GeneID[gene_clusters$Cluster %in% c("cluster4", "cluster6")], 
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       keyType = "SYMBOL",
                       pAdjustMethod ="BH",
                       readable = TRUE,
                       universe = rownames(micro_mashr$result$lfsr))
gobp_c4_c6 <- gobp_c4_c6@result



# gobp_all <- enrichGO(gene = gene_clusters$GeneID[gene_clusters$Cluster %in% c("cluster4", "cluster6", "cluster1", "cluster2")],
#                      OrgDb = org.Hs.eg.db,
#                      ont = "BP",
#                      keyType = "SYMBOL",
#                      pAdjustMethod ="BH",
#                      readable = TRUE,
#                      universe = rownames(micro_mashr$result$lfsr))
# gobp_all <- gobp_all@result

####################

saveRDS(gobp_c1_c2, "./analysis/microglia/differential_expression/microglia_gobp_c1_c2.rds")
saveRDS(gobp_c4_c6, "./analysis/microglia/differential_expression/microglia_gobp_c4_c6.rds")
saveRDS(gobp_all, "./analysis/microglia/differential_expression/microglia_gobp_all.rds")

########################################
########################################
########################################

# rank genes for GSEA

micro_mashr <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_mashr_object.rds")
all_lfsrs <- as.data.frame(micro_mashr$result$lfsr)
all_betas <- as.data.frame(micro_mashr$result$PosteriorMean)

####################

log10lfsr <- as.data.frame(-log10(all_lfsrs) * sign(all_betas))

max <- max(as.numeric(unlist(log10lfsr))[is.finite(as.numeric(unlist(log10lfsr)))])
min <- min(as.numeric(unlist(log10lfsr))[is.finite(as.numeric(unlist(log10lfsr)))])

log10lfsr <- replace(log10lfsr, log10lfsr > max, max * 10)
log10lfsr <- replace(log10lfsr, log10lfsr < min, min * 10)


CR <- log10lfsr[c("Mathys_EC_CR", "Mathys_HIP_CR", "Mathys_FC_PFC_CR", "Mathys_TC_MTG_CR", "Mathys_AnG_CR", "Mathys_TH_CR")]
CR$median <- apply(CR, 1, median)


non_CR <- log10lfsr[,!grepl("CR", colnames(log10lfsr))]
non_CR$median <- apply(non_CR, 1, median)


t <- data.frame(gene = rownames(CR),
                CR_med = CR$median,
                non_CR_med = non_CR$median)
t$diff <- t$CR_med - t$non_CR_med
t <- t %>%
  arrange(desc(diff))

########################################

list <- t$diff
names(list) <- t$gene

gse <- gseGO(geneList = list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH",
             eps = 0)

########################################

list_kegg <- t$diff
names(list_kegg) <- mapIds(org.Hs.eg.db, keys=t$gene, column="ENTREZID", keytype="SYMBOL")
list_kegg <- list_kegg[!is.na(names(list_kegg))]


kegg <- gseKEGG(geneList = list_kegg,
                keyType = "kegg", 
                minGSSize = 3, 
                maxGSSize = 800, 
                verbose = TRUE, 
                organism = "hsa", 
                pAdjustMethod = "BH",
                eps = 0)

####################

saveRDS(gse, "./analysis/microglia/differential_expression/microglia_GSEA_GO.rds")
saveRDS(kegg, "./analysis/microglia/differential_expression/microglia_GSEA_kegg.rds")

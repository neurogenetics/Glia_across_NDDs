library(tidyverse)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(corrplot)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

cts <- readRDS("./analysis/cellculture/iMG_bulk_RNA/cleaned_counts/all_batches_cts_merged.rds")
meta <- readRDS("./analysis/cellculture/iMG_bulk_RNA/metadata/sample_meta_merged.rds")

##################################################
##################################################
##################################################

# PCA w/ all samples together

##################################################

# start by filtering genes
# filtering so genes are expressed @ cpm > 1 in at least 3 samples of any treatment group

dge <- DGEList(cts)
dge <- calcNormFactors(dge, method = "TMM")
cpm <- as.data.frame(edgeR::cpm(dge))


meta <- meta[match(colnames(cpm), rownames(meta)), ]

genes_keep <- vector("logical", nrow(cpm))
genes_per_group <- data.frame(treatment = character(), n_genes = numeric(), stringsAsFactors = FALSE)

for (t in unique(meta$treatment)) {
  samples_in_group <- rownames(meta)[meta$treatment == t]
  
  cpm_sub <- cpm[, samples_in_group, drop = FALSE]
  
  gene_pass <- rowSums(cpm_sub > 1) >= 3
  
  genes_keep <- genes_keep | gene_pass
  
  genes_per_group <- rbind(genes_per_group,
                           data.frame(treatment = t, n_genes = sum(gene_pass), stringsAsFactors = FALSE))
}

gene_pass <- rownames(cpm)[gene_pass]

cts_filtered <- cts[rownames(cts) %in% gene_pass, ]
cpm_filtered <- cpm[rownames(cpm) %in% gene_pass, ]

##################################################

pca <- prcomp(t(cpm_filtered), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pca$x)

pca <- pca %>%
  dplyr::select(PC1, PC2, PC3) %>%
  rownames_to_column(var = "sample_ID")

meta <- meta %>%
  rownames_to_column(var = "sample_ID") %>%
  left_join(pca, by = "sample_ID")

meta$PPMI_line <- as.factor(meta$PPMI_line)



ggplot(meta, aes(x = PC1, y = PC2, colour = PPMI_line)) + 
  geom_point(size = 2)


ggplot(meta, aes(x = PC1, y = PC2, colour = treatment_dose)) + 
  geom_point(size = 2)


ggplot(meta, aes(x = PC1, y = PC2, colour = differentiation)) + 
  geom_point(size = 2)

##################################################
##################################################
##################################################

# PCA on batches separately

##################################################

batch_list <- c("XR10585", "DA10634")

meta_list <- list()
cts_list <- list()
cpm_list <- list()

# filter for expressed genes in each batch

for (batch in batch_list){
  meta_filtered <- meta[grepl(batch, rownames(meta)), ]
  cts_filtered <- cts[,colnames(cts) %in% rownames(meta_filtered)]

  dge <- DGEList(cts_filtered)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm <- as.data.frame(edgeR::cpm(dge))
  
  meta_filtered <- meta_filtered[match(colnames(cpm), rownames(meta_filtered)), ]
  
  genes_keep <- vector("logical", nrow(cpm))
  genes_per_group <- data.frame(treatment = character(), n_genes = numeric(), stringsAsFactors = FALSE)
  
  for (t in unique(meta_filtered$treatment)) {
    samples_in_group <- rownames(meta_filtered)[meta_filtered$treatment == t]
    
    cpm_sub <- cpm[, samples_in_group, drop = FALSE]
    
    gene_pass <- rowSums(cpm_sub > 1) >= 3
    
    genes_keep <- genes_keep | gene_pass
    
    genes_per_group <- rbind(genes_per_group,
                             data.frame(treatment = t, n_genes = sum(gene_pass), stringsAsFactors = FALSE))
  }
  
  gene_pass <- rownames(cpm)[gene_pass]
  
  cts_filtered <- cts_filtered[rownames(cts_filtered) %in% gene_pass, ]
  cpm_filtered <- cpm[rownames(cpm) %in% gene_pass, ]
  
  meta_list[[batch]] <- meta_filtered
  cts_list[[batch]] <- cts_filtered
  cpm_list[[batch]] <- cpm_filtered
}

##################################################

# run PCA and correlate to metadata variables

PCA_corr_list <- list()
PCA_variance_explained <- list()

for (batch in batch_list){
  cpm <- cpm_list[[batch]]
  meta <- meta_list[[batch]]
  
  pca_obj <- prcomp(t(cpm), center = TRUE, scale. = TRUE)
  pca <- as.data.frame(pca_obj$x)
  pca <- pca %>%
    rownames_to_column(var = "sample")
  
  # get % variance explained
  var_explained <- summary(pca_obj)$importance[2, ] * 100 
  var_df <- tibble(
    PC = names(var_explained),
    variance_explained = var_explained
  )
  PCA_variance_explained[[batch]] <- var_df  # Save it
  
  meta <- meta %>%
    rownames_to_column(var = "sample") %>%
    left_join(pca[,1:21], by = "sample")
  
  # correlate each PC w/ vars of interest
  pc_cols <- grep("^PC\\d+$", names(meta), value = TRUE)
  vars_to_test <- c("PPMI_line", "cells_plated", "treatment", "treatment_dose_scaled", "sex", "RNA_nanodrop",
                    "pct.intronic", "pct.ribo", "pct.mito")
  
  results <- list()
  
  for (pc in pc_cols) {
    pc_results <- list()
    
    for (var in vars_to_test) {
      x <- meta[[pc]]
      y_raw <- meta[[var]]
      
      # convert categorical to numeric encoding
      if (is.factor(y_raw) || is.character(y_raw)) {
        y <- as.numeric(factor(y_raw))
      } else {
        y <- y_raw
      }
      
      # compute Pearson correlation
      res <- cor.test(x, y, method = "pearson")
      
      pc_results[[var]] <- list(
        estimate = res$estimate,
        p.value = res$p.value
      )
    }
    
    results[[pc]] <- pc_results
  }
  
  # make results df
  tidy_results <- map_dfr(names(results), function(pc) {
    map_dfr(results[[pc]], function(r) {
      tibble(
        estimate = r$estimate,
        p.value = r$p.value
      )
    }, .id = "variable") %>%
      mutate(PC = pc)
  })
  
  PCA_corr_list[[batch]] <- tidy_results
}

##################################################

# make corr heatmaps

XR10585_PC_corr <- PCA_corr_list$XR10585
DA10634_PC_corr <- PCA_corr_list$DA10634

#########################

XR10585_PC_corr_mat <- XR10585_PC_corr %>%
  dplyr::select(-p.value) %>%
  pivot_wider(names_from = PC, values_from = estimate) %>%
  column_to_rownames(var = "variable") %>%
  as.matrix()

XR10585_PC_corr_mat <- XR10585_PC_corr_mat*XR10585_PC_corr_mat


col_fun <- colorRamp2(breaks = c(0, 1), colors = c("white", "red"))
Heatmap(XR10585_PC_corr_mat,
        cluster_rows = F, 
        cluster_columns = F,
        col = col_fun,
        name = "r^2")

#########################

DA10634_PC_corr_mat <- DA10634_PC_corr %>%
  dplyr::select(-p.value) %>%
  pivot_wider(names_from = PC, values_from = estimate) %>%
  column_to_rownames(var = "variable") %>%
  as.matrix()

DA10634_PC_corr_mat <- DA10634_PC_corr_mat*DA10634_PC_corr_mat


Heatmap(DA10634_PC_corr_mat,
        cluster_rows = F, 
        cluster_columns = F,
        col = col_fun,
        name = "r^2")

##################################################
##################################################
##################################################

# plot % variance explained

XR10585_var <- PCA_variance_explained$XR10585 %>%
  arrange(desc(variance_explained))
XR10585_var$PC <- factor(XR10585_var$PC, levels = XR10585_var$PC)

ggplot(XR10585_var, aes(x = PC, y = variance_explained)) + 
  geom_point()



DA10634_var <- PCA_variance_explained$DA10634 %>%
  arrange(desc(variance_explained))
DA10634_var$PC <- factor(DA10634_var$PC, levels = DA10634_var$PC)

ggplot(DA10634_var, aes(x = PC, y = variance_explained)) + 
  geom_point()

##################################################
##################################################
##################################################

# plotting % variance explained for a PC against correlation of treatment w/ that PC

t <- XR10585_PC_corr %>%
  subset(variable == "treatment_dose_scaled") %>%
  left_join(XR10585_var, by = "PC") %>%
  arrange(desc(variance_explained)) %>%
  mutate(r2 = estimate*estimate) %>%
  subset(p.value > 0.1 & variance_explained > 1)

t$PC <- factor(t$PC, levels = t$PC)

ggplot(t, aes(x = r2, y = variance_explained)) + 
  geom_point(aes(colour = PC))

ggplot(t, aes(x = variable, y = r2)) + 
  geom_point() +
  facet_wrap(~PC) +
  theme(axis.text.x = element_text(angle = 90))

#########################

t <- XR10585_PC_corr %>%
  left_join(XR10585_var, by = "PC") %>%
  arrange(desc(variance_explained)) %>%
  mutate(r2 = estimate*estimate)
  
t2 <- t %>%
  subset(variable %in% c("treatment", "treatment_dose_scaled")) %>%
  group_by(PC) %>%
  summarise(sum_tx = sum(r2))

t3 <- t %>%
  left_join(t2, by = "PC") %>%
  subset(!variable %in% c("treatment", "treatment_dose_scaled")) %>%
  group_by(PC) %>%
  mutate(max_other = max(r2)) %>%
  mutate(diff = max_other - sum_tx) %>%
  filter(diff > 0) %>%
  pull(PC) %>%
  unique()

t4 <- t %>%
  subset(variable %in% c("treatment", "treatment_dose_scaled")) %>%
  subset(p.value > 0.1 & variance_explained > 1) %>%
  group_by(PC) %>%
  summarise(n = n()) %>%
  filter(n == 2) %>%
  pull(PC)

PC_final <- t4[!t4 %in% setdiff(t4, t3)]






t$PC <- factor(t$PC, levels = t$PC)

ggplot(t, aes(x = r2, y = variance_explained)) + 
  geom_point(aes(colour = PC))

ggplot(t, aes(x = variable, y = r2)) + 
  geom_point() +
  facet_wrap(~PC) +
  theme(axis.text.x = element_text(angle = 90))



##################################################
##################################################
##################################################

# checking what covars are co-linear w/ each other

meta_XR10585 <- meta %>%
  rownames_to_column(var = "sample") %>%
  mutate(treatment_dose = as.numeric(factor(treatment_dose))) %>%
  subset(grepl("XR10585", sample))

corrplot(cor(meta_XR10585[, c(3, 10, 12:14, 16)]), col = colorRampPalette(c("#05409e", "white", "#e64a02"))(200))



meta_DA10634 <- meta %>%
  rownames_to_column(var = "sample") %>%
  mutate(treatment_dose = as.numeric(factor(treatment_dose))) %>%
  subset(grepl("DA10634", sample))

corrplot(cor(meta_DA10634[, c(3, 10, 12:14, 16)]), col = colorRampPalette(c("#05409e", "white", "#e64a02"))(200))
















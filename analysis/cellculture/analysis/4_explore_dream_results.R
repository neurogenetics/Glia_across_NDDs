library(tidyverse)
library(fgsea)
library(eulerr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(ComplexHeatmap)
library(circlize)
library(purrr)

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
DA10682 <- readRDS("./analysis/cellculture/iMG_bulk_RNA/tmp_results/DA10682_dream_res_df_list.rds")
DA10725 <- readRDS("./analysis/cellculture/iMG_bulk_RNA/tmp_results/DA10725_dream_res_df_list.rds")



dream_res_list <- c(XR10585, DA10634,DA10682,DA10725)

##################################################
##################################################
##################################################

# count number of significant DEGs

nsigs_df <- data.frame()

for (test in names(dream_res_list)){
  df <- dream_res_list[[test]]
  
  n_sigs <- as.numeric(nrow(df[df$adj.P.Val < 0.05 & abs(df$logFC) > 0.5, ]))
  
  nsigs_df <- rbind(nsigs_df, 
                    data.frame(test = test,
                               n_sigs = n_sigs))
}

##################################################
##################################################
##################################################

# filter all dfs for DAM genes + make heatmap

DAM_res <- list()
logFC_list <- list()

for (test in names(dream_res_list)){
  df <- dream_res_list[[test]]
  
  DAM_df <- df[rownames(df) %in% DAM_genes, ]
  
  DAM_res[[test]] <- DAM_df
  
  
  df1 <- DAM_df %>%
    rownames_to_column(var = "gene") %>%
    dplyr::select(gene, logFC) %>%
    dplyr::rename(!!paste0(test) := logFC)
  
  logFC_list[[test]] <- df1
}

DAM_res_df <- purrr::reduce(logFC_list, full_join, by = "gene")


DAM_res_df <- DAM_res_df %>%
  column_to_rownames(var = "gene") %>%
  # dplyr::select(-`Vorinostat_1e-05`, -`Entinostat_1e-08`) %>%
  as.matrix()

DAM_res_df <- DAM_res_df[, 
                         # grepl("UV-DAMPS|Entinostat|Vorinostat", colnames(DAM_res_df)) #for figure on previously demonstrated
                         # grepl("IFNg|Spe|Ibr|Vac|Api", colnames(DAM_res_df)) #for figure on novel
                         ]%>%na.omit()

Heatmap(DAM_res_df,
        cluster_columns = F)

##################################################
##################################################
##################################################

# comparing high dose vorinostat to lower dose

high <- dream_res_list$`Entinostat_1e-06`
med <- dream_res_list$`Entinostat_1e-07`
low <- dream_res_list$`Entinostat_1e-08`

high <- high %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(logFC_high = logFC)

low <- low %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(logFC_low = logFC)

summary <- med %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(logFC_med = logFC) %>%
  left_join(high, by = "gene") %>%
  left_join(low, by = "gene")

ggplot(summary, aes(x = logFC_low, y = logFC_med)) +
  geom_point() + 
  geom_abline()

ggplot(summary, aes(x = logFC_med, y = logFC_high)) +
  geom_point() + 
  geom_abline()

##################################################

summary_filtered <- summary[summary$gene %in% DAM_genes, ]

##################################################
##################################################
##################################################

# comparing our vorinostat to PDJ vorinostat
# for iMG, they used 0.1uM = 1e-7M, our low dose
# they present "blue" differentiation in the paper

PDJ_blue <- read.csv("./analysis/iMicroglia/test/iMG_blue.csv")
PDJ_orange <- read.csv("./analysis/iMicroglia/test/iMG_orange.csv")
PDJ_HMC3_vor <- read.csv("./analysis/iMicroglia/test/HMC3_vorinostat.csv")

PDJ_blue <- PDJ_blue %>%
  dplyr::rename(logFC_PDJ_blue = log2FoldChange, gene = GeneID) %>%
  dplyr::select(gene, logFC_PDJ_blue)

PDJ_orange <- PDJ_orange %>%
  dplyr::rename(logFC_PDJ_orange = log2FoldChange, gene = GeneID) %>%
  dplyr::select(gene, logFC_PDJ_orange)

PDJ_HMC3_vor <- PDJ_HMC3_vor %>%
  dplyr::rename(logFC_PDJ_HMC3 = log2FoldChange, gene = GeneID) %>%
  dplyr::select(gene, logFC_PDJ_HMC3)

summary <- summary %>%
  left_join(PDJ_blue, by = "gene") %>%
  left_join(PDJ_orange, by = "gene") %>%
  left_join(PDJ_HMC3_vor, by = "gene") %>%
  na.omit()

ggplot(summary, aes(x = logFC_med, y = logFC_PDJ_blue)) +
  geom_point() + 
  geom_abline()

##################################################

# filter df for any condition which has logfc >|1|

genes_PDJ <- read.csv("./analysis/iMicroglia/test/genes_PDJ_supptable.csv", header = F)
genes_PDJ <- genes_PDJ$V1

# summary_filtered <- summary[rowSums(abs(summary[,-1]) > 1) >= 1, ]
summary_filtered <- summary[summary$gene %in% DAM_genes, ]
# summary_filtered <- summary[summary$gene %in% genes_PDJ, ]

beta_vals <- summary_filtered[, -1]  # Remove gene names
gene_names <- summary_filtered[, 1]

# Get all unique pairs of columns
combs <- combn(colnames(beta_vals), 2, simplify = FALSE)

# Initialize result holder
results <- data.frame(
  comparison = character(),
  same_sign = integer(),
  discordant = integer(),
  zero_involved = integer(),
  stringsAsFactors = FALSE
)

# Compare each pair
for (pair in combs) {
  col1 <- beta_vals[[pair[1]]]
  col2 <- beta_vals[[pair[2]]]
  
  # Sign comparison
  sign1 <- sign(col1)
  sign2 <- sign(col2)
  
  # Zero handling
  zero_involved <- sum(sign1 == 0 | sign2 == 0)
  
  # Same sign (+1/+1 or -1/-1)
  same_sign <- sum((sign1 == sign2) & (sign1 != 0))
  
  # Discordant (+1/-1 or -1/+1)
  discordant <- sum((sign1 * sign2) == -1)
  
  # Store result
  results <- rbind(results, data.frame(
    comparison = paste(pair[1], "vs", pair[2]),
    same_sign = same_sign,
    discordant = discordant,
    zero_involved = zero_involved
  ))
}

# sort(genes_PDJ[!genes_PDJ %in% summary_filtered$gene])


p1 <- summary_filtered$gene[sign(summary_filtered$logFC_med) == sign(summary_filtered$logFC_low)]
p2 <- summary_filtered$gene[sign(summary_filtered$logFC_med) == sign(summary_filtered$logFC_PDJ_blue)]
p3 <- summary_filtered$gene[sign(summary_filtered$logFC_med) == sign(summary_filtered$logFC_high)]

summary(p1 %in% p2)
intersect(p2, p3)
p1[!p1 %in% p3]
sort(p1[p1 %in% p2])

##################################################








##################################################
##################################################
##################################################




low <- rownames(dream_res$IFNg_2.5)[dream_res$`IFNg_2.5`$adj.P.Val <= 0.05 & abs(dream_res$`IFNg_2.5`$logFC) > 0]
med <- rownames(dream_res$`IFNg_12.5`)[dream_res$`IFNg_12.5`$adj.P.Val <= 0.05 & abs(dream_res$`IFNg_12.5`$logFC) > 0]
high <- rownames(dream_res$`IFNg_25`)[dream_res$`IFNg_25`$adj.P.Val <= 0.05 & abs(dream_res$`IFNg_25`$logFC) > 0]

plot(venn(list("low" = low, "med" = med, "high" = high)))

x = union(low, union(med, high))



low_df <- dream_res$`IFNg_2.5`[rownames(dream_res$`IFNg_2.5`) %in% x, ]
med_df <- dream_res$`IFNg_12.5`[rownames(dream_res$`IFNg_12.5`) %in% x, ]
high_df <- dream_res$`IFNg_25`[rownames(dream_res$`IFNg_25`) %in% x, ]

low_df <- low_df %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(`low` = logFC)

med_df <- med_df %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(`med` = logFC)

high_df <- high_df %>%
  rownames_to_column(var = "gene") %>%
  dplyr::select(gene, logFC) %>%
  dplyr::rename(`high` = logFC) %>%
  left_join(low_df, by = "gene") %>%
  left_join(med_df, by = "gene")

df <- high_df %>%
  pivot_longer(cols = c(`high`, `med`, `low`), names_to = "treatment")

ggplot(df, aes(x = value, colour = treatment)) +
  geom_density()

##################################################
##################################################
##################################################

# apoptosis: hsa04210

pathway_id <- "hsa04210"

# Retrieve pathway data
pathway_data <- keggGet(pathway_id)

# Extract Entrez Gene IDs
genes_raw <- pathway_data[[1]]$GENE
entrez_ids <- genes_raw[seq(1, length(genes_raw), 2)]

# Convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")



fgsea_res <- data.frame()

for (test in names(dream_res_list)){
  test_df <- dream_res_list[[test]]
  
  list <- test_df$z.std
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


list_test <- dream_res_list$`UV-DAMPS_336000`$z.std
names(list_test) <- rownames(dream_res_list$`UV-DAMPS_336000`)

plotEnrichment(pathway = DAM_genes,
               stats = list_test)


# leaving off 070825 -- add IFNg sig + make dotplot

##################################################

EnhancedVolcano(dream_res$`PFF_1e-06`$logFC,
                lab = rownames(dream_res$`PFF_1e-06`$logFC),
                x = 'logFC',
                y = 'adj.P.Val')


ggplot(dream_res$`PFF_1e-06`, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.58) +
  geom_vline(xintercept = -0.58)



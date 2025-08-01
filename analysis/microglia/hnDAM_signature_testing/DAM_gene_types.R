library(tidyverse)
library(rtracklayer)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################

# get gene types for DAM genes

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")


gtf <- import("/data/ADRD/human_brain_atlasing/3_fastq_processing/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz")
gtf_df <- as.data.frame(gtf)

gtf_df <- gtf_df[!duplicated(gtf_df$gene_name), ]
gtf_df <- gtf_df[gtf_df$gene_name %in% DAM_genes, ]

summary(factor(gtf_df$gene_type))

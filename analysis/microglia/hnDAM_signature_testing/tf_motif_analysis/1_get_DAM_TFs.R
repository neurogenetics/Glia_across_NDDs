library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

hocomoco_v11 <- read.csv("./analysis/microglia/DAM_signature_discovery/tf_motif_analysis/HUMAN_mono_motifs.tsv", sep = "\t")
hocomoco_v11$Transcription.factor <- gsub("HUMAN:", "", hocomoco_v11$Transcription.factor)

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")

t <- DAM_genes[DAM_genes %in% hocomoco_v11$Transcription.factor]

# MITF, PPARG, ARID5B, MAFB
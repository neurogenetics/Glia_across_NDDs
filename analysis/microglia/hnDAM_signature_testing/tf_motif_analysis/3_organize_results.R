library(tidyverse)
library(fgsea)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

gene_map <- read.csv("./analysis/microglia/DAM_signature_discovery/tf_motif_analysis/cellranger_features.tsv", sep = "\t", header = F)

gene_map <- gene_map %>%
  dplyr::rename(ensembl_id = V1, gene = V2) %>%
  dplyr::select(-V3)

##################################################
##################################################
##################################################

MAFB <- read.csv("./analysis/microglia/DAM_signature_discovery/tf_motif_analysis/outs/fimo_MAFB.H13CORE.0.PS.A_meme_format/MAFB.H13CORE.0.PS.A_meme_format_motif_hits.tsv",
                 header = F, sep = "\t")
MITF <- read.csv("./analysis/microglia/DAM_signature_discovery/tf_motif_analysis/outs/fimo_MITF.H13CORE.0.P.B_meme_format/MITF.H13CORE.0.P.B_meme_format_motif_hits.tsv",
                 header = F, sep = "\t")
PPARG <- read.csv("./analysis/microglia/DAM_signature_discovery/tf_motif_analysis/outs/fimo_PPARG.H13CORE.0.P.B_meme_format/PPARG.H13CORE.0.P.B_meme_format_motif_hits.tsv",
                  header = F, sep = "\t")

##################################################

# convert ensembl id to gene name

MAFB <- MAFB %>%
  dplyr::rename(ensembl_id = V10) %>%
  left_join(gene_map, by = "ensembl_id")

MITF <- MITF %>%
  dplyr::rename(ensembl_id = V10) %>%
  left_join(gene_map, by = "ensembl_id")

PPARG <- PPARG %>%
  dplyr::rename(ensembl_id = V10) %>%
  left_join(gene_map, by = "ensembl_id")

##################################################

# group & sort by number of motifs found per gene

MAFB_hits <- MAFB %>%
  group_by(gene) %>%
  summarise(n_motifs = n()) %>%
  arrange(desc(n_motifs)) %>%
  rownames_to_column(var = "rank")

MITF_hits <- MITF %>%
  group_by(gene) %>%
  summarise(n_motifs = n()) %>%
  arrange(desc(n_motifs)) %>%
  rownames_to_column(var = "rank")

PPARG_hits <- PPARG %>%
  group_by(gene) %>%
  summarise(n_motifs = n()) %>%
  arrange(desc(n_motifs)) %>%
  rownames_to_column(var = "rank")

##################################################

# get gene lists by n_motifs

MAFB_genes <- MAFB_hits$gene[1:300]

MITF_genes <- MITF_hits$gene[1:300]
# MITF_genes <- MITF_hits$gene[MITF_hits$n_motifs > 50]

PPARG_genes <- PPARG_hits$gene[1:300]

pathways <- list(MAFB = MAFB_genes,
                 MITF = MITF_genes,
                 PPARG = PPARG_genes)

##################################################
##################################################
##################################################

DAM_genes_ranking <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_GEPs_zscore_df.rds")
DAM_genes_rank <- DAM_genes_ranking$avg_z
names(DAM_genes_rank) <- DAM_genes_ranking$gene


plotEnrichment(pathway = pathways$MITF, stats = DAM_genes_rank)


res = fgsea(pathways = pathways, stats = DAM_genes_rank)


intersect(PPARG_genes, DAM_genes_ranking$gene[1:105])


t <- MITF_hits[MITF_hits$gene %in% DAM_genes_ranking$gene[DAM_genes_ranking$avg_z > 1], ]







library(tidyverse)
library(fgsea)
library(DESeq2)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
pathway <- list(DAM_genes)
names(pathway) <- "DAM"

##################################################
##################################################
##################################################

# Ruzicka 2024, Science (schizophrenia study, not able to access data so using DE results from supp table 4)

res_SCZ_Ruzicka <- read.csv("./analysis/microglia/DAM_signature_replication/non_NDDs/Ruzicka_2024/Ruzicka_2024_micro_DEGs.csv")

#########################

res_SCZ_Ruzicka$t <- sign(res_SCZ_Ruzicka$MtSinai_logFC) * -log10(res_SCZ_Ruzicka$MtSinai_p_val)
res_SCZ_Ruzicka <- res_SCZ_Ruzicka %>%
  arrange(desc(t))

ranking_SCZ_Ruzicka <- res_SCZ_Ruzicka$t
names(ranking_SCZ_Ruzicka) <- res_SCZ_Ruzicka$gene

gsea_SCZ_Ruzicka <- fgsea(pathways = pathway, stats = ranking_SCZ_Ruzicka, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_SCZ_Ruzicka)

##################################################
##################################################
##################################################

# Frohlich 2024, Nat. Neuro (aging & neuropsychiatric disorder study, DEGs from supp. tables 6 and 16, respectively)

##################################################

# aging

res_aging_Frohlich <- read.csv("./analysis/microglia/DAM_signature_replication/non_NDDs/Frohlich_2024/Frohlich_2024_aging_micro_DEGs.csv")

res_aging_Frohlich$t <- sign(res_aging_Frohlich$logFC) * -log10(res_aging_Frohlich$P.Value)
res_aging_Frohlich <- res_aging_Frohlich %>%
  arrange(desc(t))

ranking_aging_Frohlich <- res_aging_Frohlich$t
names(ranking_aging_Frohlich) <- res_aging_Frohlich$ID

gsea_aging_Frohlich <- fgsea(pathways = pathway, stats = ranking_aging_Frohlich, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_aging_Frohlich)

##################################################

# NPS

res_NPS_Frohlich <- read.csv("./analysis/microglia/DAM_signature_replication/non_NDDs/Frohlich_2024/Frohlich_2024_NPS_micro_DEGs.csv")

res_NPS_Frohlich$t <- sign(res_NPS_Frohlich$logFC) * -log10(res_NPS_Frohlich$P.Value)
res_NPS_Frohlich <- res_NPS_Frohlich %>%
  arrange(desc(t))

ranking_NPS_Frohlich <- res_NPS_Frohlich$t
names(ranking_NPS_Frohlich) <- res_NPS_Frohlich$ID

gsea_NPS_Frohlich <- fgsea(pathways = pathway, stats = ranking_NPS_Frohlich, scoreType = "std")
plotEnrichment(pathway = pathway$DAM, stats = ranking_NPS_Frohlich)

library(tidyverse)
library(mashr)

setwd("/data/ADRD/glia_across_NDDs")

########################################

micro_props_glm_res <- readRDS("./analysis/microglia/cluster_proportions/micro_props_logit_GLM_results_1v1.rds")
micro_props_glm_res$term <- gsub("group", "", micro_props_glm_res$term)
micro_props_glm_res <- micro_props_glm_res %>%
  mutate(region = sub(".*_", "", d_r))
micro_props_glm_res$disease_region <- paste0(micro_props_glm_res$term, "_", micro_props_glm_res$region)

########################################

# create shat and bhat for SE and estimates

bhat <- micro_props_glm_res %>%
  pivot_wider(names_from = "disease_region", id_cols = "cluster", values_from = "Estimate") %>%
  column_to_rownames(var = "cluster") %>%
  as.matrix()

shat <- micro_props_glm_res %>%
  pivot_wider(names_from = "disease_region", id_cols = "cluster", values_from = "Std. Error") %>%
  column_to_rownames(var = "cluster") %>%
  as.matrix()

########################################

# create mashed dataset
mash_data = mash_set_data(bhat, shat)

########################################

# set up & mash canonical cov matrix
U.c = cov_canonical(mash_data)

m.c = mash(mash_data, U.c)

########################################

# create & mash data-driven cov matrices
m.1by1 = mash_1by1(mash_data)

strong = get_significant_results(m.1by1, thresh = 0.05)

U.pca = cov_pca(mash_data, npc = 4, subset = strong)

U.ed = cov_ed(mash_data, U.pca, subset = strong)

m.ed = mash(mash_data, c(U.c, U.ed))

saveRDS(m.ed, "./analysis/microglia/cluster_proportions/micro_props_mashr_obj_1v1.rds")

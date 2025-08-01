library(tidyverse)
library(mashr)

setwd("/data/ADRD/glia_across_NDDs")

########################################

oligo_props_glm_res <- readRDS("./analysis/oligodendrocytes/cluster_proportions/oligo_props_logit_GLM_results_1v1.rds")
oligo_props_glm_res$term <- gsub("group", "", oligo_props_glm_res$term)
oligo_props_glm_res <- oligo_props_glm_res %>%
  mutate(region = sub(".*_", "", d_r))
oligo_props_glm_res$disease_region <- paste0(oligo_props_glm_res$term, "_", oligo_props_glm_res$region)

########################################

# create shat and bhat for SE and estimates

bhat <- oligo_props_glm_res %>%
  pivot_wider(names_from = "disease_region", id_cols = "cluster", values_from = "Estimate") %>%
  column_to_rownames(var = "cluster") %>%
  as.matrix()

shat <- oligo_props_glm_res %>%
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

saveRDS(m.c, "./analysis/oligodendrocytes/cluster_proportions/oligo_props_mashr_obj_canonical.rds")

library(tidyverse)
library(mashr)

setwd("/data/ADRD/glia_across_NDDs")

########################################

astro_props_glm_res <- readRDS("./analysis/astrocytes/cluster_proportions/cortical_astro_props_logit_GLM_results_1v1.rds")
astro_props_glm_res$term <- gsub("group", "", astro_props_glm_res$term)
astro_props_glm_res <- astro_props_glm_res %>%
  mutate(region = sub(".*_", "", d_r))
astro_props_glm_res$disease_region <- paste0(astro_props_glm_res$term, "_", astro_props_glm_res$region)

########################################

# create shat and bhat for SE and estimates

bhat <- astro_props_glm_res %>%
  pivot_wider(names_from = "disease_region", id_cols = "cluster", values_from = "Estimate") %>%
  column_to_rownames(var = "cluster") %>%
  as.matrix()

shat <- astro_props_glm_res %>%
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

saveRDS(m.c, "./analysis/astrocytes/cluster_proportions/cortical_astros_props_mashr_obj_canonical.rds")

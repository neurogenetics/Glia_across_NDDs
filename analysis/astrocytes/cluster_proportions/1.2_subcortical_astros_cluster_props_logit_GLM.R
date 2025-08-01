library(tidyverse)
library(car)
library(lme4)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

astro_celllevel_meta <- readRDS("./analysis/astrocytes/subcortical_astros_metadata_ANNOTATED.rds")
donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")

##################################################

# cell counts & props by donor x region

astro_celllevel_meta$donor_region <- paste0(astro_celllevel_meta$donor, "_", astro_celllevel_meta$region)
t <- as.data.frame(table(astro_celllevel_meta$cluster_anno, astro_celllevel_meta$donor_region))

props <- t %>%
  group_by(Var2) %>%
  mutate(total_donor_cells = sum(Freq), cluster_proportion = Freq / total_donor_cells) %>%
  ungroup() %>%
  rename(donor_region = Var2) %>%
  rename(cluster = Var1) %>%
  rename(cell_counts_in_cluster = Freq) %>%
  extract(donor_region, into = c("donor", "region"), regex = "(.+)_(.+)") %>%
  left_join(donor_meta, by = "donor")

saveRDS(props, "./analysis/astrocytes/cluster_proportions/subcortical_astro_raw_donor_region_props.rds")

########################################
########################################
########################################

# GLM analysis of astroglia proportions w/ logit transformation

########################################

astro_props <- readRDS("./analysis/astrocytes/cluster_proportions/subcortical_astro_raw_donor_region_props.rds")
donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")


astro_props$dataset_region <- paste0(astro_props$study, "_", astro_props$region)
dataset_regions <- unique(astro_props$dataset_region)

# dataset_region <- dataset_regions[1]

########################################

overall_GLM_summary <- data.frame()

for (dataset_region in dataset_regions){
  t <- astro_props[astro_props$dataset_region == dataset_region, ]
  
  diseases <- unique(t$group[!t$group %in% "HC"])
  
  for (disease in diseases){
    t2 <- t[t$group %in% c("HC", disease), ]
    
    mat <- t2 %>%
      dplyr::select(donor, cluster_proportion, cluster) %>%
      pivot_wider(names_from = donor, values_from = cluster_proportion) %>%
      column_to_rownames(var = "cluster")
    
    epsilon <- 1e-6
    mat[mat == 0] <- epsilon
    mat[mat == 1] <- 1 - epsilon
    
    mat_trans <- logit(mat)
    
    prop_long <- mat_trans %>%
      rownames_to_column("cluster") %>%
      pivot_longer(-cluster, names_to = "donor", values_to = "prop") %>%
      left_join(donor_meta, by = "donor") %>%
      mutate(group = factor(group, levels = c("HC", setdiff(unique(group), "HC"))))
    
    prop_long$age_scaled <- scale(prop_long$age)
    
    glm_results <- prop_long %>%
      group_by(cluster) %>%
      summarise(
        model = list(glm(prop ~ group + sex + age, data = cur_data(), family = gaussian())),
        .groups = "drop")
    
    glm_summary <- glm_results %>%
      mutate(coef_summary = map(model, ~summary(.x)$coefficients),
             p_values = map(coef_summary, ~as.data.frame(.x) %>% 
                              rownames_to_column("term"))) %>%
      dplyr::select(cluster, p_values) %>%
      unnest(p_values) %>% 
      filter(str_detect(term, "group")) 
    
    glm_summary$d_r <- dataset_region
    
    overall_GLM_summary <- rbind(overall_GLM_summary, glm_summary)
  }
}

astro_glm_summary_combined <- overall_GLM_summary %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  mutate(p_adj = p.adjust(`Pr(>|t|)`, method = "fdr"))

astro_glm_summary_combined <- astro_glm_summary_combined[!astro_glm_summary_combined$`Pr(>|t|)` == "NaN", ]

saveRDS(astro_glm_summary_combined, file = "./analysis/astrocytes/cluster_proportions/subcortical_astro_props_logit_GLM_results_1v1.rds")

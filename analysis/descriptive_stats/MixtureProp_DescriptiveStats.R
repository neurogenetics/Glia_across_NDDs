library(tidyverse)
library(mashr)

setwd("/data/ADRD/glia_across_NDDs")
micro_mashr <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_mashr_object.rds")
astro_mashr <- readRDS("./analysis/astrocytes//differential_expression/astrocytes_diffexp_mashr_object.rds")
oligo_mashr <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_mashr_object.rds")

#Micro: MI contributions
pl = data.frame(mi = get_estimated_pi(micro_mashr))
pl$group = rownames(pl)
pl <- pl %>%
  arrange(desc(mi))
sum(pl$mi[grepl("^(Gerrits|Pineda|Mathys|NM|AMP)", pl$group)])#need to inclue AMP and NM since not renamed in this simple load-in

#ASTRO: MI contributions
pl = data.frame(mi = get_estimated_pi(astro_mashr))
pl$group = rownames(pl)
pl <- pl %>%
  arrange(desc(mi))
sum(pl$mi[grepl("^(Gerrits|Pineda|Mathys|NM|AMP)", pl$group)])#need to inclue AMP and NM since not renamed in this simple load-in

#Oligo: # MI contributions
pl = data.frame(mi = get_estimated_pi(oligo_mashr))
pl$group = rownames(pl)
pl <- pl %>%
  arrange(desc(mi))
sum(pl$mi[grepl("^(Gerrits|Pineda|Mathys|NM|AMP)", pl$group)])#need to inclue AMP and NM since not renamed in this simple load-in


#####ASTRO: now for specific comparisons from PD|ALS how many comparisons have corr > 0.3 with FTD|FTLD
tPCA_corr <- astro_mashr$fitted_g$Ulist$ED_tPCA
#Get variable names that contain "PD" or "ALS"
pd_als_vars <- grep("PD|ALS", rownames(tPCA_corr), value = TRUE)
ftd_vars <- grep("FTD|FTLD", colnames(tPCA_corr), value = TRUE)
subset_matrix <- tPCA_corr[pd_als_vars, ftd_vars, drop = FALSE]
# Step 3: Find PD/ALS variables that have any covariance > 0.3 with FTD variables
vars_with_high_covar <- rownames(subset_matrix)[apply(subset_matrix, 1, function(x) any(x > 0.3))]
count_vars <- length(vars_with_high_covar)


#####OLIGO: now for specific comparisons from ALS how many comparisons have corr > 0.3 with !ALS
tPCA_corr <- oligo_mashr$fitted_g$Ulist$ED_tPCA
als_vars <- grep("ALS", rownames(tPCA_corr), value = TRUE)
non_als_vars <- grep("ALS", colnames(tPCA_corr), value = TRUE, invert = TRUE)
# Subset the matrix to ALS rows and non-ALS columns
subset_matrix <- tPCA_corr[als_vars, non_als_vars, drop = FALSE]
vars_with_high_covar <- colnames(subset_matrix)[apply(subset_matrix, 1, function(x) any(x > 0.3))]
length(vars_with_high_covar)

#####MICROS: now for specific comparisons from ALS how many comparisons have corr > 0.3 with !ALS
# Load correlation matrix
tPCA_corr <- micro_mashr$fitted_g$Ulist$ED_tPCA
als_vars <- grep("ALS", rownames(tPCA_corr), value = TRUE)
non_als_vars <- grep("ALS", colnames(tPCA_corr), value = TRUE, invert = TRUE)
# Subset the matrix to ALS rows and non-ALS columns
subset_matrix <- tPCA_corr[als_vars, non_als_vars, drop = FALSE]
vars_with_high_covar <- colnames(subset_matrix)[apply(subset_matrix, 1, function(x) any(x > 0.3))]
length(vars_with_high_covar)

#####Astro: now for specific comparisons from ALS how many comparisons have corr > 0.3 with !ALS
# Load correlation matrix
tPCA_corr <- astro_mashr$fitted_g$Ulist$ED_tPCA
als_vars <- grep("ALS", rownames(tPCA_corr), value = TRUE)
non_als_vars <- grep("ALS", colnames(tPCA_corr), value = TRUE, invert = TRUE)
# Subset the matrix to ALS rows and non-ALS columns
subset_matrix <- tPCA_corr[als_vars, non_als_vars, drop = FALSE]
vars_with_high_covar <- colnames(subset_matrix)[apply(subset_matrix, 1, function(x) any(x > 0.3))]
length(vars_with_high_covar)


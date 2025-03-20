library(tidyverse)

PD_meta <- readRDS("/data/ADRD/glia_across_NDDs/amp_pd/combined_metadata_prefilter.rds")
ALS_meta <- readRDS("/data/ADRD/glia_across_NDDs/alsftd/combined_metadata_prefilter.rds")
AD_meta <- readRDS("/data/ADRD/glia_across_NDDs/ad/combined_metadata_prefilter.rds")
FTD_meta <- readRDS("/data/ADRD/glia_across_NDDs/gerrits/combined_metadata_prefilter.rds")

####################

# PD

round(mean(PD_meta$nCount_RNA)) # 11377
t1 <- as.data.frame(quantile(PD_meta$nCount_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t1 <- rownames_to_column(t1, var = "percentile")


round(mean(PD_meta$nFeature_RNA)) # 3696
t2 <- as.data.frame(quantile(PD_meta$nFeature_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t2 <- rownames_to_column(t2, var = "percentile")

mean(PD_meta$pct.mito) # 0.501806
mean(PD_meta$pct.ribo) # 0.6415221

####################

# ALS

round(mean(ALS_meta$nCount_RNA)) # 8775
t3 <- as.data.frame(quantile(ALS_meta$nCount_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t3 <- rownames_to_column(t3, var = "percentile")

round(mean(ALS_meta$nFeature_RNA)) # 2967
t4 <- as.data.frame(quantile(ALS_meta$nFeature_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t4 <- rownames_to_column(t4, var = "percentile")

mean(ALS_meta$pct.mito) # 2.071297
mean(ALS_meta$pct.ribo) # 0.8006054

####################

# AD

round(mean(AD_meta$nCount_RNA)) # 10807
t5 <- as.data.frame(quantile(AD_meta$nCount_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t5 <- rownames_to_column(t5, var = "percentile")

round(mean(AD_meta$nFeature_RNA)) # 3390
t6 <- as.data.frame(quantile(AD_meta$nFeature_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t6 <- rownames_to_column(t6, var = "percentile")

mean(AD_meta$pct.mito) # 2.961997
mean(AD_meta$pct.ribo) # 0.6798433

####################

# Gerrits FTD

t13 <- as.data.frame(quantile(FTD_meta$nCount_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t13 <- rownames_to_column(t13, var = "percentile")

t14 <- as.data.frame(quantile(FTD_meta$nFeature_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t14 <- rownames_to_column(t14, var = "percentile")

####################

counts_all <- left_join(t1, t3, by = "percentile") %>%
  left_join(t5, by = "percentile")
colnames(counts_all) <- c("percentile", "PD", "ALS", "AD")


features_all <- left_join(t2, t4, by = "percentile") %>%
  left_join(t6, by = "percentile")
colnames(features_all) <- c("percentile", "PD", "ALS", "AD")

####################

# set filters & examine data post-filter

PD_meta_filtered <- PD_meta %>%
  filter(nCount_RNA < 43277 & nCount_RNA > 2062 & nFeature_RNA > 1333 & pct.mito < 2 & pct.ribo < 2)
t7 <- as.data.frame(quantile(PD_meta_filtered$nCount_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t7 <- rownames_to_column(t7, var = "percentile")
t8 <- as.data.frame(quantile(PD_meta_filtered$nFeature_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t8 <- rownames_to_column(t8, var = "percentile")


ALS_meta_filtered <- ALS_meta %>%
  filter(nCount_RNA < 41758 & nCount_RNA > 812 & nFeature_RNA > 619 & pct.mito < 2 & pct.ribo < 2)
t9 <- as.data.frame(quantile(ALS_meta_filtered$nCount_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t9 <- rownames_to_column(t9, var = "percentile")
t10 <- as.data.frame(quantile(ALS_meta_filtered$nFeature_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t10 <- rownames_to_column(t10, var = "percentile")


AD_meta_filtered <- AD_meta %>%
  filter(nCount_RNA < 46886 & nCount_RNA > 825 & nFeature_RNA > 581 & pct.mito < 2 & pct.ribo < 2)
t11 <- as.data.frame(quantile(AD_meta_filtered$nCount_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t11 <- rownames_to_column(t11, var = "percentile")
t12 <- as.data.frame(quantile(AD_meta_filtered$nFeature_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1)))
t12 <- rownames_to_column(t12, var = "percentile")



FTD_meta_filtered <- FTD_meta %>%
  filter(nCount_RNA < 8182 & nCount_RNA > 465 & nFeature_RNA > 383.0 & pct.mito < 2 & pct.ribo < 2)


counts_all_postfilter <- left_join(t7, t9, by = "percentile") %>%
  left_join(t11, by = "percentile")
colnames(counts_all_postfilter) <- c("percentile", "PD", "ALS", "AD")


features_all_postfilter <- left_join(t8, t10, by = "percentile") %>%
  left_join(t12, by = "percentile")
colnames(features_all_postfilter) <- c("percentile", "PD", "ALS", "AD")


####################


saveRDS(PD_meta_filtered, file = "/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/PD_metadata_FILTERED.rds")
saveRDS(ALS_meta_filtered, file = "/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/ALSFTD_metadata_FILTERED.rds")
saveRDS(AD_meta_filtered, file = "/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/AD_metadata_FILTERED.rds")
saveRDS(FTD_meta_filtered, file = "/data/ADRD/glia_across_NDDs/combined_data/post_filter_metadata/FTD_metadata_FILTERED.rds")

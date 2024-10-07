library(Seurat)
library(tidyverse)
library(scCustomize)

fastq_list <- read.csv(file = "/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt", sep = "\t", header = F)
fastq_list <- fastq_list$V1

# loop through all the cellranger_seurat RDS files, pull relevant metadata

combined_metadata <- data.frame()

for (fastq in fastq_list){
  cellranger_seurat_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/", fastq, "/", fastq, "_cellranger_seurat.rds")
  cellranger_seurat <- readRDS(cellranger_seurat_dir)
  
  meta <- cellranger_seurat@meta.data[, c("donor", "batch", "region", "nCount_RNA", "nFeature_RNA", "pct.mito", "pct.ribo")]
  
  combined_metadata <- rbind(combined_metadata, meta)
}

# filter out cells from samples that had the same donor twice in the same pool

t <- as.data.frame(table(combined_metadata$batch))

combined_metadata_filtered <- combined_metadata %>%
  filter(!(batch == "PM-PD_Set4_1" & donor == "PM-MS_70978")) %>%
  filter(!(batch == "PM-PD_Set4_2" & donor == "PM-MS_70978")) %>%
  filter(!(batch == "PM-PD_Set5_1" & donor == "PM-MS_46475")) %>%
  filter(!(batch == "PM-PD_Set5_2" & donor == "PM-MS_46475")) %>%
  filter(!(batch == "PM-PD_Set7a_1" & donor == "PM-MS_55245")) %>%
  filter(!(batch == "PM-PD_Set7a_2" & donor == "PM-MS_55245"))

t2 <- as.data.frame(table(combined_metadata_filtered$batch))

saveRDS(combined_metadata_filtered, file = "/data/ADRD/glia_across_NDDs/amp_pd/combined_metadata_prefilter.rds")

########################################################################

meta <- readRDS("/data/ADRD/glia_across_NDDs/amp_pd/outs/combined_metadata_prefilter.rds")

# get QC metric stats

# gene count
nCount_mean <- round(mean(meta$nCount_RNA)) # 11377
quantile(meta$nCount_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1))
# 2.5%        5%       10%       25%       50%       75%       90%       95%     97.5%      100% 
# 2062.00   2642.00   3394.00   4994.00   7758.00  13170.00  23867.00  33489.00  43277.05 574882.00 

# number of features
nFeat_mean <- round(mean(meta$nFeature_RNA)) # 3696
quantile(meta$nFeature_RNA, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1))
# 2.5%      5%     10%     25%     50%     75%     90%     95%   97.5%    100% 
# 1333.0  1596.9  1904.0  2463.0  3248.0  4466.0  6299.0  7479.0  8361.0 19103.0 

# pct mito
mito_mean <- mean(meta$pct.mito) # 0.5018

# pct ribo
ribo_mean <- mean(meta$pct.ribo) # 0.6415


# get number of cells per donor by batch
meta$batch <- gsub("PM-PD_", "", meta$batch)
meta$donor_batch <- paste0(meta$donor, "_", meta$batch)
t3 <- as.data.frame(table(meta$donor_batch))
t3 <- t3 %>% 
  arrange(Freq) %>% 
  select(Var1, Freq)

ggplot(t3, aes(x = 1:nrow(t3), y = Freq)) + 
  geom_point(size = 0.1) +
  ggtitle("# cells per donor x pool") +
  labs(x = "donor",
       y = "# cells")

# filter donor x batch w/ less than 400 cells
donors_to_keep <- as.character(t3[54:nrow(t3), "Var1"])


# filter
meta_filtered <- meta %>%
  filter(nCount_RNA < 43277 & nFeature_RNA > 1333 & pct.mito < 2 & pct.ribo < 2) %>%
  filter(donor_batch %in% donors_to_keep)

######

saveRDS(meta_filtered, file = "/data/ADRD/glia_across_NDDs/amp_pd/combined_metadata_postfilter.rds")

########################################################################

meta_filtered <- readRDS("/data/ADRD/glia_across_NDDs/amp_pd/combined_metadata_postfilter.rds")
cells_to_keep <- rownames(meta_filtered)

# filter CellBender seurat objects to only have cells to keep

meta_filtered$fastq <- paste0("PM-PD_", meta_filtered$batch)
fastq_list_filtered <- unique(meta_filtered$fastq)

for (fastq in fastq_list_filtered){
  cellbender_seurat_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/", fastq, "/", fastq, "_cellbender_seurat.rds")
  cellbender_seurat <- readRDS(cellbender_seurat_dir)
  
  cells_in_sample <- Cells(cellbender_seurat)
  cells_to_keep_in_sample <- intersect(cells_in_sample, cells_to_keep)
  
  cellbender_seurat_filtered <- subset(cellbender_seurat, cells = cells_to_keep_in_sample)
  
  filtered_save_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/", fastq, "/", fastq, "_cellbender_seurat_filtered.rds")
  saveRDS(cellbender_seurat_filtered, file = filtered_save_dir)
}

########################################################################

meta_filtered <- readRDS("/data/ADRD/glia_across_NDDs/amp_pd/combined_metadata_postfilter.rds")
meta_filtered$fastq <- paste0("PM-PD_", meta_filtered$batch)
fastq_list_filtered <- unique(meta_filtered$fastq)

# load each filtered object, subset by the regions, then save region objects into lists to merge

DMV_list <- list()
GPi_list <- list()
M1_list <- list()
DLPFC_list <- list()
V1_list <- list()

regions <- c("DMNX", "GPI", "PMC", "PFC", "PVC")

# fastq_list_test <- c("PM-PD_Set78_C1", "PM-PD_Set60_E2")

for (fastq in fastq_list_filtered){
  cellbender_filtered_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/", fastq, "/", fastq, "_cellbender_seurat_filtered.rds")
  cellbender_seurat_filtered <- readRDS(cellbender_filtered_dir)
  
  regions_in_obj <- unique(cellbender_seurat_filtered$region)
  
  for (brainregion in regions){
    if (brainregion %in% regions_in_obj) {
      fastq_brainregion <- subset(cellbender_seurat_filtered, subset = region == brainregion)
      
      if (brainregion == "DMNX") {
        DMV_list[[length(DMV_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "GPI") {
        GPi_list[[length(GPi_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "PMC") {
        M1_list[[length(M1_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "PFC") {
        DLPFC_list[[length(DLPFC_list) + 1]] <- fastq_brainregion
      } else if (brainregion == "PVC") {
        V1_list[[length(V1_list) + 1]] <- fastq_brainregion
      }
    }
  }
}


# merge regions into their own objects

DMV_merged <- Merge_Seurat_List(DMV_list)
saveRDS(DMV_merged, file = "/data/ADRD/glia_across_NDDs/amp_pd/DMV_merged.rds")
rm(DMV_merged)

GPi_merged <- Merge_Seurat_List(GPi_list)
saveRDS(GPi_merged, file = "/data/ADRD/glia_across_NDDs/amp_pd/GPi_merged.rds")
rm(GPi_merged)

M1_merged <- Merge_Seurat_List(M1_list)
saveRDS(M1_merged, file = "/data/ADRD/glia_across_NDDs/amp_pd/M1_merged.rds")
rm(M1_merged)

V1_merged <- Merge_Seurat_List(V1_list)
saveRDS(V1_merged, file = "/data/ADRD/glia_across_NDDs/amp_pd/V1_merged.rds")
rm(V1_merged)

DLPFC_merged <- Merge_Seurat_List(DLPFC_list)
saveRDS(DLPFC_merged, file = "/data/ADRD/glia_across_NDDs/amp_pd/DLPFC_merged.rds")
rm(DLPFC_merged)

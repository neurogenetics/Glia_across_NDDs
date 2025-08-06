library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

files <- list.files("./analysis/cellculture/iMG_bulk_RNA/fastq/batch1_XR10585/", full.names = F)

R1 <- files[grepl("_R1_", files)]
R2 <- files[grepl("_R2_", files)]

samplesheet <- data.frame(sample = NA,
                          fastq_1 = paste0("/data/ADRD/glia_across_NDDs/analysis/cellculture/iMG_bulk_RNA/fastq/batch1_XR10585/", R1[1]),
                          fastq_2 = paste0("/data/ADRD/glia_across_NDDs/analysis/cellculture/iMG_bulk_RNA/fastq/batch1_XR10585/", R2[1]),
                          strandedness = "auto")

samplesheet$sample <- sub("_.*", "", R1[1])

write.csv(samplesheet, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/1_indexing_run/samplesheet.csv", row.names = FALSE)

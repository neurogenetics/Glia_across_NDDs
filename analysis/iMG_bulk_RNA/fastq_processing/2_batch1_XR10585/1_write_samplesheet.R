library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

files <- list.files("./analysis/cellculture/iMG_bulk_RNA/fastq/batch1_XR10585/", full.names = F)

R1 <- files[grepl("_R1_", files)]
R2 <- files[grepl("_R2_", files)]

samplesheet <- data.frame(sample = NA,
                          fastq_1 = paste0("/data/ADRD/glia_across_NDDs/analysis/cellculture/iMG_bulk_RNA/fastq/batch1_XR10585/", R1),
                          fastq_2 = paste0("/data/ADRD/glia_across_NDDs/analysis/cellculture/iMG_bulk_RNA/fastq/batch1_XR10585/", R2),
                          strandedness = "auto")

samplesheet$sample <- sub("_.*", "", R1)


samplesheet_1 <- samplesheet[1:54,]
samplesheet_2 <- samplesheet[55:108,]


write.csv(samplesheet_1, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/2_batch1_XR10585/samplesheet_1.csv", row.names = FALSE)
write.csv(samplesheet_2, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/2_batch1_XR10585/samplesheet_2.csv", row.names = FALSE)

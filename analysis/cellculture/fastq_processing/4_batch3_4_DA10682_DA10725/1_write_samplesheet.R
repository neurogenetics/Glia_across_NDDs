library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

files <- list.files("./analysis/cellculture/iMG_bulk_RNA/fastq/batch3_4_DA10682_DA10725/", full.names = F)

R1 <- files[grepl("_R1_", files)]
R2 <- files[grepl("_R2_", files)]

samplesheet <- data.frame(sample = NA,
                          fastq_1 = paste0("/data/ADRD/glia_across_NDDs/analysis/cellculture/iMG_bulk_RNA/fastq/batch3_4_DA10682_DA10725/", R1),
                          fastq_2 = paste0("/data/ADRD/glia_across_NDDs/analysis/cellculture/iMG_bulk_RNA/fastq/batch3_4_DA10682_DA10725/", R2),
                          strandedness = "auto")

samplesheet$sample <- sub("_.*", "", R1)


samplesheet_1 <- samplesheet[1:96,]
samplesheet_2 <- samplesheet[97:192,]
samplesheet_3 <- samplesheet[192:294,] ## THIS WAS A MISTAKE!!! should have been 193:294, need to remove sample 141 from batch 3 processing
samplesheet_4 <- samplesheet[295:396,]
samplesheet_5 <- samplesheet[397:498,]
samplesheet_6 <- samplesheet[499:600,]
samplesheet_7 <- samplesheet[601:702,]
samplesheet_8 <- samplesheet[703:804,]
samplesheet_9 <- samplesheet[805:888,]
samplesheet_10 <- samplesheet[889:941,]


write.csv(samplesheet_1, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_1.csv", row.names = FALSE)
write.csv(samplesheet_2, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_2.csv", row.names = FALSE)
write.csv(samplesheet_3, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_3.csv", row.names = FALSE)
write.csv(samplesheet_4, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_4.csv", row.names = FALSE)
write.csv(samplesheet_5, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_5.csv", row.names = FALSE)
write.csv(samplesheet_6, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_6.csv", row.names = FALSE)
write.csv(samplesheet_7, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_7.csv", row.names = FALSE)
write.csv(samplesheet_8, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_8.csv", row.names = FALSE)
write.csv(samplesheet_9, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_9.csv", row.names = FALSE)
write.csv(samplesheet_10, "./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_10.csv", row.names = FALSE)

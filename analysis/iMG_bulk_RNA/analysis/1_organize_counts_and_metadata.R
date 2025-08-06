library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

##################################################
##################################################
##################################################

# read in counts data

cts_1 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch1_XR10585/processing_batch1/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_2 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch1_XR10585/processing_batch2/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_3 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch2_DA10634/processing_batch1/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_4 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch2_DA10634/processing_batch2/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")

#########################

cts_1$gene_name <- make.unique(cts_1$gene_name)
cts_1 <- cts_1 %>%
  dplyr::select(-gene_id) %>%
  column_to_rownames(var = "gene_name")

cts_2$gene_name <- make.unique(cts_2$gene_name)
cts_2 <- cts_2 %>%
  dplyr::select(-gene_id) %>%
  column_to_rownames(var = "gene_name")

cts_3$gene_name <- make.unique(cts_3$gene_name)
cts_3 <- cts_3 %>%
  dplyr::select(-gene_id) %>%
  column_to_rownames(var = "gene_name")

cts_4$gene_name <- make.unique(cts_4$gene_name)
cts_4 <- cts_4 %>%
  dplyr::select(-gene_id) %>%
  column_to_rownames(var = "gene_name")


identical(rownames(cts_1), rownames(cts_2))
identical(rownames(cts_2), rownames(cts_3))
identical(rownames(cts_3), rownames(cts_4))


cts_merged <- cbind(cts_1, cts_2, cts_3, cts_4)

cts_merged <- round(cts_merged)
colnames(cts_merged) <- gsub("\\.", "_", colnames(cts_merged))

#########################

saveRDS(cts_merged, file = "./analysis/cellculture/iMG_bulk_RNA/cleaned_counts/all_batches_cts_merged.rds")

##################################################
##################################################
##################################################

# organize metadata and add technical stuff from multiqc

##################################################

cts <- readRDS("./analysis/cellculture/iMG_bulk_RNA/cleaned_counts/all_batches_cts_merged.rds")

#########################

sequencing_meta <- read.csv("./analysis/cellculture/iMG_bulk_RNA/metadata/iMG_RNA_metadata.csv")

##################################################

# get qualimap data (pct intronic)

qualimap_1 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch1_XR10585/processing_batch1/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_2 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch1_XR10585/processing_batch2/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_3 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch2_DA10634/processing_batch1/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_4 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch2_DA10634/processing_batch2/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")



qualimap <- rbind(qualimap_1, qualimap_2, qualimap_3, qualimap_4)

qualimap$pct.intronic <- qualimap$Intronic / rowSums(qualimap[2:4]) * 100
qualimap <- qualimap %>%
  dplyr::select(Sample, pct.intronic) %>%
  dplyr::rename(sample_ID = Sample)

qualimap$sample_ID <- gsub("-", "_", qualimap$sample_ID)

##################################################

# add sex info for lines

sequencing_meta$sex <- ifelse(sequencing_meta$PPMI_line %in% c("3480", "3453"), "F", "M")

##################################################

# calculate pct.mito / pct.ribo

mito_genes <- grep("^MT-", rownames(cts), value = TRUE, ignore.case = TRUE)
ribo_genes <- grep("^RP[SL]", rownames(cts), value = TRUE, ignore.case = TRUE)

total_counts <- colSums(cts)

pct.mito <- (colSums(cts[mito_genes, , drop = FALSE]) / total_counts) * 100
pct.ribo <- (colSums(cts[ribo_genes, , drop = FALSE]) / total_counts) * 100

df <- data.frame(sample_ID = colnames(cts),
                 pct.mito = pct.mito,
                 pct.ribo = pct.ribo)

rownames(df) <- NULL

##################################################

meta_merged <- sequencing_meta %>%
  left_join(qualimap, by = "sample_ID") %>%
  left_join(df, by = "sample_ID")

meta_merged$sample_ID <- gsub("-", "_", meta_merged$sample_ID)

meta_merged <- meta_merged %>%
  column_to_rownames(var = "sample_ID")

meta_merged <- meta_merged[colnames(cts), ]

meta_merged$treatment_dose_scaled <- paste0(meta_merged$treatment, "_", meta_merged$dose_scaled)
meta_merged$treatment_dose <- paste0(meta_merged$treatment, "_", meta_merged$dose)

meta_merged <- na.omit(meta_merged)

#########################

saveRDS(meta_merged, file = "./analysis/cellculture/iMG_bulk_RNA/metadata/sample_meta_merged.rds")

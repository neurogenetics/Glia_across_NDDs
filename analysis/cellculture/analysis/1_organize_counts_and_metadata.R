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
cts_5 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch1/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_6 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch2/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_7 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch3/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_8 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch4/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_9 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch5/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_10 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch6/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_11 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch7/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_12 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch8/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_13 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch9/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
cts_14 <- read.csv("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch10/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")

# accidentally included sample 141 in two batches--should stay in cts_6, remove from cts_7
cts_7 <- cts_7[, !grepl("141", colnames(cts_7))]


cts_list <- list(cts_1, cts_2, cts_3, cts_4, cts_5, 
                 cts_6, cts_7, cts_8, cts_9, cts_10,
                 cts_11, cts_12, cts_13, cts_14)

#########################

cts_list <- lapply(cts_list, function(cts) {
  cts$gene_name <- make.unique(cts$gene_name)
  cts %>%
    select(-gene_id) %>%
    column_to_rownames(var = "gene_name")
})


cts_merged <- do.call(cbind, cts_list)

cts_merged <- round(cts_merged)
colnames(cts_merged) <- gsub("DA.10682", "DA10682", colnames(cts_merged))
colnames(cts_merged) <- gsub("DA.10725", "DA10725", colnames(cts_merged))
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

sequencing_meta$sample_ID <- gsub("DA-10634", "DA10634", sequencing_meta$sample_ID)
sequencing_meta$sample_ID <- gsub("DA-10725", "DA10725", sequencing_meta$sample_ID)
sequencing_meta$sample_ID <- gsub("-", "_", sequencing_meta$sample_ID)

##################################################

# get qualimap data (pct intronic)

qualimap_1 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch1_XR10585/processing_batch1/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_2 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch1_XR10585/processing_batch2/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_3 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch2_DA10634/processing_batch1/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_4 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch2_DA10634/processing_batch2/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_5 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch1/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_6 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch2/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_7 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch3/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_8 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch4/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_9 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch5/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_10 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch6/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_11 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch7/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_12 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch8/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_13 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch9/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")
qualimap_14 <- read.delim("./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch10/multiqc/star_salmon/multiqc_report_data/qualimap_genomic_origin.txt")


# accidentally included sample 141 in two batches--should stay in qualimap_6, remove from qualimap_7
qualimap_7 <- qualimap_7[!grepl("141", qualimap_7$Sample),]


qualimap <- rbind(qualimap_1, qualimap_2, qualimap_3, qualimap_4, qualimap_5,
                  qualimap_6, qualimap_7, qualimap_8, qualimap_9, qualimap_10,
                  qualimap_11, qualimap_12, qualimap_13, qualimap_14)
                  
qualimap$pct.intronic <- qualimap$Intronic / rowSums(qualimap[2:4]) * 100
qualimap$total_reads <- rowSums(qualimap[2:4])
qualimap <- qualimap %>%
  dplyr::select(Sample, pct.intronic, total_reads) %>%
  dplyr::rename(sample_ID = Sample)


qualimap$sample_ID <- gsub("DA-10682", "DA10682", qualimap$sample_ID)
qualimap$sample_ID <- gsub("DA-10725", "DA10725", qualimap$sample_ID)
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

meta_merged$treatment_dose <- paste0(meta_merged$treatment, "_", meta_merged$dose)

# meta_merged <- na.omit(meta_merged)

#########################

saveRDS(meta_merged, file = "./analysis/cellculture/iMG_bulk_RNA/metadata/sample_meta_merged.rds")

library(tidyverse)
library(scCustomize)
library(Seurat)

args <- commandArgs(trailingOnly = T)
fastq <- args[1]

# (for testing) fastq <- "PM-PD_Set83_E1"

####################################

demuxafy_combined_results_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/", fastq, "/combined/combined_results.tsv")
combined_results <- read.delim(demuxafy_combined_results_dir)


# filter for cells that both Dropulation and Demuxalot call singlets (remove ALL heterogenic doublets)
filtered <- combined_results[combined_results$Vireo_DropletType == "singlet" &
                               combined_results$Demuxalot_DropletType == "singlet", ]


# filter out cells that Dropulation and Demuxalot assign to different donors
filtered <- filtered[filtered$Vireo_Individual_Assignment == 
                       filtered$Demuxalot_Individual_Assignment, ]


# retain cells which the majority (2/3) of the general doublet detectors call singlets
filtered <- filtered[rowSums(filtered[, c("DoubletDetection_DropletType", 
                                          "scDblFinder_DropletType", 
                                          "scds_DropletType")] == "singlet") >= 2, ]


# filter for droplets which CellBender calls real
cellbender_barcodes_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellbender/outs/", fastq, "/cellbender_cell_barcodes.csv")
cellbender_barcodes <- read.csv(cellbender_barcodes_dir, header = F)

filtered <- filtered[filtered$Barcode %in% cellbender_barcodes$V1, ]

####################################

# read in raw CellRanger output and make Seurat object
cellranger_mat_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/", fastq, "/raw_feature_bc_matrix/")
cellranger_mat <- Read10X(cellranger_mat_dir)
cellranger_seurat <- CreateSeuratObject(counts = cellranger_mat)


# filter for droplets which pass QC for demux, doublets, and ambient RNA
cellranger_seurat <- cellranger_seurat[, colnames(cellranger_seurat) %in% filtered$Barcode]


# UMI/percent MT/ribo filtering
cellranger_seurat[["pct.mito"]] <- PercentageFeatureSet(cellranger_seurat, pattern = "^MT-")
cellranger_seurat[["pct.ribo"]] <- PercentageFeatureSet(cellranger_seurat, pattern = "^RP")


# add relevant metadata
meta <- filtered[filtered$Barcode %in% colnames(cellranger_seurat), ]
meta <- meta[c("Barcode", "Demuxalot_Individual_Assignment")]
colnames(meta) <- c("barcode", "donor")
batch <- gsub("^.*_(Set\\d+_C\\d+)$", "\\1", fastq)
meta$barcode_batch <- paste0(meta$barcode, "_", batch)
meta$batch <- batch

cellranger_meta <- cellranger_seurat@meta.data
cellranger_meta <- rownames_to_column(cellranger_meta, var = "barcode")


# add brain regions
regions_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/", fastq, "/donor_region_list.txt")
regions <- read.delim(regions_dir, header = F, sep = "\t")

split <- do.call(rbind, strsplit(regions$V1, "(?<=.)-(?!.*-)", perl = TRUE))
regions <- data.frame(donor = split[, 1], region = split[, 2])

meta <- meta %>%
  left_join(regions, by = "donor")


# combine all metadata
combined_cellranger_meta <- cellranger_meta %>%
  left_join(meta, by = "barcode") %>%
  column_to_rownames(var = "barcode")

identical(cellranger_meta$barcode, rownames(combined_cellranger_meta))

cellranger_seurat@meta.data <- combined_cellranger_meta

####################################

# do the same for the CellBender corrected H5
cellbender_mat_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellbender/outs/", fastq, "/cellbender_filtered.h5")
cellbender_mat <- Read_CellBender_h5_Mat(cellbender_mat_dir)
cellbender_seurat <- CreateSeuratObject(counts = cellbender_mat, names.field = 1, names.delim = "_")

cellbender_seurat <- cellbender_seurat[, colnames(cellbender_seurat) %in% colnames(cellranger_seurat)]


cellbender_meta <- cellbender_seurat@meta.data
cellbender_meta <- rownames_to_column(cellbender_meta, var = "barcode")

combined_cellbender_meta <- cellbender_meta %>%
  left_join(meta, by = "barcode") %>%
  column_to_rownames(var = "barcode")

identical(cellbender_meta$barcode, rownames(combined_cellbender_meta))

cellbender_seurat@meta.data <- combined_cellbender_meta

####################################

# make sure that both objects have exactly the same cells (already filtered the other way above)
cellranger_seurat <- cellranger_seurat[, colnames(cellranger_seurat) %in% colnames(cellbender_seurat)]

# make the barcodes appended w/ the batch
colnames(cellranger_seurat) <- cellranger_seurat$barcode_batch
colnames(cellbender_seurat) <- cellbender_seurat$barcode_batch

####################################

cellranger_seurat_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/", fastq, "/", fastq, "_cellranger_seurat.rds")
cellbender_seurat_dir <- paste0("/data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/", fastq, "/", fastq, "_cellbender_seurat.rds")


# check that cell IDs match between objects and save objects if they match
if (!all(colnames(cellranger_seurat) %in% colnames(cellbender_seurat)) || 
    !all(colnames(cellbender_seurat) %in% colnames(cellranger_seurat))) {
  stop("Error: The cell IDs do not match between cellranger_seurat and cellbender_seurat.")
} else {
  print("All cell IDs match!")
  saveRDS(cellranger_seurat, file = cellranger_seurat_dir)
  saveRDS(cellbender_seurat, file = cellbender_seurat_dir)
}

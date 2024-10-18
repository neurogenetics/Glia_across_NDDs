library(tidyverse)
library(scCustomize)
library(Seurat)

args <- commandArgs(trailingOnly = T)
fastq <- args[1]

# (for testing) fastq <- "SRR27882078"

####################################

demuxafy_combined_results_dir <- paste0("/data/ADRD/ALSFTD_multiregion/fastq_processing/demuxafy/outs/", fastq, "/combined/combined_results.tsv")
combined_results <- read.delim(demuxafy_combined_results_dir)


# retain cells which the majority (2/3) of the general doublet detectors call singlets
filtered <- combined_results[rowSums(combined_results[, c("DoubletDetection_DropletType",
                                                          "scDblFinder_DropletType",
                                                          "scds_DropletType")] == "singlet") >= 2, ]


# filter for droplets which CellBender calls real
cellbender_barcodes_dir <- paste0("/data/ADRD/ALSFTD_multiregion/fastq_processing/cellbender/outs/", fastq, "/cellbender_cell_barcodes.csv")
cellbender_barcodes <- read.csv(cellbender_barcodes_dir, header = F)

filtered <- filtered[filtered$Barcode %in% cellbender_barcodes$V1, ]

####################################

# read in raw CellRanger output and make Seurat object
cellranger_mat_dir <- paste0("/data/ADRD/ALSFTD_multiregion/fastq_processing/cellranger/outs/", fastq, "/raw_feature_bc_matrix/")
cellranger_mat <- Read10X(cellranger_mat_dir)
cellranger_seurat <- CreateSeuratObject(counts = cellranger_mat)


# filter for droplets which pass QC for demux, doublets, and ambient RNA
cellranger_seurat <- cellranger_seurat[, colnames(cellranger_seurat) %in% filtered$Barcode]


# UMI/percent MT/ribo filtering
cellranger_seurat[["pct.mito"]] <- PercentageFeatureSet(cellranger_seurat, pattern = "^MT-")
cellranger_seurat[["pct.ribo"]] <- PercentageFeatureSet(cellranger_seurat, pattern = "^RP")

####################################

# do the same for the CellBender corrected H5
cellbender_mat_dir <- paste0("/data/ADRD/ALSFTD_multiregion/fastq_processing/cellbender/outs/", fastq, "/cellbender_filtered.h5")
cellbender_mat <- Read_CellBender_h5_Mat(cellbender_mat_dir)
cellbender_seurat <- CreateSeuratObject(counts = cellbender_mat, names.field = 1, names.delim = "_")

cellbender_seurat <- cellbender_seurat[, colnames(cellbender_seurat) %in% colnames(cellranger_seurat)]

####################################

# make sure that both objects have exactly the same cells (already filtered the other way above)
cellranger_seurat <- cellranger_seurat[, colnames(cellranger_seurat) %in% colnames(cellbender_seurat)]

####################################

# for converting b/w fastq name and sample name
meta <- read.csv(file = "/data/ADRD/ALSFTD_multiregion/fastq_processing/basic_meta.csv", header = T)
meta <- meta[meta$fastq %in% fastq, ]
region <- meta$region
sample <- substring(meta$sample, 1, 3)


cellranger_seurat_dir <- paste0("/data/ADRD/ALSFTD_multiregion/fastq_processing/final_outs/", fastq, "/", sample, "_", region, "_cellranger_seurat.rds")
cellbender_seurat_dir <- paste0("/data/ADRD/ALSFTD_multiregion/fastq_processing/final_outs/", fastq, "/", sample, "_", region, "_cellbender_seurat.rds")


# check that cell IDs match between objects and save objects if they match
if (!all(colnames(cellranger_seurat) %in% colnames(cellbender_seurat)) || 
    !all(colnames(cellbender_seurat) %in% colnames(cellranger_seurat)) ||
    ncol(cellranger_seurat) > ncol(cellbender_seurat) ||
    ncol(cellranger_seurat) < ncol(cellbender_seurat)
    )  {
  stop("Error: The cell IDs do not match between cellranger_seurat and cellbender_seurat.")
} else {
  print("All cell IDs match!")
  saveRDS(cellranger_seurat, file = cellranger_seurat_dir)
  saveRDS(cellbender_seurat, file = cellbender_seurat_dir)
}


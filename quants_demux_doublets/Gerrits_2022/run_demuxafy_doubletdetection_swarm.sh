#!/bin/bash

FASTQ=${1}
COUNTS=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix
BARCODES=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix/barcodes.tsv.gz
DOUBLETDETECTION_OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/doubletdetection

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/.singularity

singularity exec --bind /data/ADRD/gerrits_ftd_snRNA/fastq_processing Demuxafy.sif DoubletDetection.py \
    -m $COUNTS \
    -b $BARCODES \
    -o $DOUBLETDETECTION_OUTDIR 

#!/bin/bash

FASTQ=${1}
COUNTS=/data/ADRD/AD_multiregion/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix
BARCODES=/data/ADRD/AD_multiregion/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix/barcodes.tsv
DOUBLETDETECTION_OUTDIR=/data/ADRD/AD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/doubletdetection

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/AD_multiregion/fastq_processing/demuxafy/.singularity

singularity exec --bind /data/ADRD/AD_multiregion/fastq_processing Demuxafy.sif DoubletDetection.py \
    -m $COUNTS \
    -b $BARCODES \
    -o $DOUBLETDETECTION_OUTDIR 

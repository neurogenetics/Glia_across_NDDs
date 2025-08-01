#!/bin/bash

FASTQ=${1}
COUNTS=/data/ADRD/ALSFTD_multiregion/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix
SCDBLFINDER_OUTDIR=/data/ADRD/ALSFTD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/scdblfinder

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/ALSFTD_multiregion/fastq_processing/demuxafy/.singularity

singularity exec --bind /data/ADRD/ALSFTD_multiregion/fastq_processing Demuxafy.sif scDblFinder.R \
    -o $SCDBLFINDER_OUTDIR \
    -t $COUNTS

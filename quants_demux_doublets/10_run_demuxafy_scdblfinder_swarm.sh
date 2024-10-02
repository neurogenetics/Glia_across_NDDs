#!/bin/bash

FASTQ=${1}
COUNTS=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix
SCDBLFINDER_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/scdblfinder

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/.singularity

singularity exec --bind /data/ADRD/amp_pd/transcriptomics/fastq_processing Demuxafy.sif scDblFinder.R \
    -o $SCDBLFINDER_OUTDIR \
    -t $COUNTS

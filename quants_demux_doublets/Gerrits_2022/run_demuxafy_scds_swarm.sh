#!/bin/bash

FASTQ=${1}
COUNTS=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix
SCDS_OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/scds

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/.singularity

singularity exec --bind /data/ADRD/gerrits_ftd_snRNA/fastq_processing Demuxafy.sif scds.R \
    -o $SCDS_OUTDIR \
    -t $COUNTS

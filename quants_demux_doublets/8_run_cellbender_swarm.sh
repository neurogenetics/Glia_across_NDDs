#!/bin/bash

FASTQ=${1}

cd /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellbender/outs/$FASTQ

INPUT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/raw_feature_bc_matrix.h5
OUTPUT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellbender/outs/$FASTQ/cellbender.h5

module load cellbender

cellbender remove-background \
    --cuda \
    --input $INPUT \
    --output $OUTPUT \
    --epochs 150

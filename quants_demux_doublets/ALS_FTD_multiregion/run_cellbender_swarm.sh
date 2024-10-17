#!/bin/bash

FASTQ=${1}

cd /data/ADRD/ALSFTD_multiregion/fastq_processing/cellbender/outs/$FASTQ

INPUT=/data/ADRD/ALSFTD_multiregion/fastq_processing/cellranger/outs/$FASTQ/raw_feature_bc_matrix.h5
OUTPUT=/data/ADRD/ALSFTD_multiregion/fastq_processing/cellbender/outs/$FASTQ/cellbender.h5

module load cellbender

cellbender remove-background \
    --cuda \
    --input $INPUT \
    --output $OUTPUT \
    --epochs 150

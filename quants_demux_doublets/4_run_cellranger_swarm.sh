#!/bin/bash

FASTQ=${1}
REFERENCE=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/inputs/refdata-gex-GRCh38-2024-A
OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ
FASTQ_DIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/$FASTQ

module load cellranger/8.0.1

cellranger count --id $FASTQ \
    --fastqs $FASTQ_DIR \
    --sample $FASTQ \
    --transcriptome $REFERENCE \
    --output-dir $OUTDIR \
    --create-bam=true \
    --jobmode=local \
    --localmem=120 \
    --localcores=$SLURM_CPUS_PER_TASK --maxjobs=20 

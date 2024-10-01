#!/bin/bash

FASTQ=${1}
REFERENCE=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/inputs/refdata-gex-GRCh38-2024-A
OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ
FASTQ_DIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/$FASTQ
TMPDIR=/lscratch/$SLURM_JOB_ID

module load cellranger/8.0.1

cellranger count --id $FASTQ \
    --fastqs $FASTQ_DIR \
    --sample $FASTQ \
    --transcriptome $REFERENCE \
    --output-dir $TMPDIR \
    --create-bam=true \
    --jobmode=local \
    --localmem=120 \
    --localcores=$SLURM_CPUS_PER_TASK --maxjobs=20 
cp $TMPDIR/outs/raw_feature_bc_matrix.h5 $OUTDIR
cp -r $TMPDIR/outs/raw_feature_bc_matrix $OUTDIR
cp $TMPDIR/outs/filtered_feature_bc_matrix.h5 $OUTDIR
cp -r $TMPDIR/outs/filtered_feature_bc_matrix $OUTDIR
cp $TMPDIR/outs/possorted_genome_bam.bam $OUTDIR
cp $TMPDIR/outs/web_summary.html $OUTDIR

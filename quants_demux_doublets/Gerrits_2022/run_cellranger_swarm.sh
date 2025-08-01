#!/bin/bash

# in order to run all samples at the same time, need to write outputs to temp dir then copy necessary outputs from lscratch to local 

FASTQ=${1}
REFERENCE=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/inputs/refdata-gex-GRCh38-2024-A
OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/cellranger/outs/$FASTQ

# make a comma-separated list of directories w/ all fastqs for each sample
FASTQ_LIST=$(find /data/ADRD/gerrits_ftd_snRNA/fastq_processing/fastqs/$FASTQ -mindepth 1 -maxdepth 1 -type d -printf "%p," | sed 's/,$//')
TMPDIR=/lscratch/$SLURM_JOB_ID

module load cellranger/8.0.1

cellranger count --id $FASTQ \
    --sample bamtofastq \
    --fastqs $FASTQ_LIST \
    --transcriptome $REFERENCE \
    --output-dir $TMPDIR \
    --create-bam=false \
    --jobmode=local \
    --localmem=120 \
    --localcores=$SLURM_CPUS_PER_TASK --maxjobs=20 
cp $TMPDIR/outs/raw_feature_bc_matrix.h5 $OUTDIR
cp -r $TMPDIR/outs/raw_feature_bc_matrix $OUTDIR
cp $TMPDIR/outs/filtered_feature_bc_matrix.h5 $OUTDIR
cp -r $TMPDIR/outs/filtered_feature_bc_matrix $OUTDIR
cp $TMPDIR/outs/web_summary.html $OUTDIR

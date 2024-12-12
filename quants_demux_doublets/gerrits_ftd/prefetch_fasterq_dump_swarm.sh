#! /bin/bash

SRR=${1}
PREFETCH_OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/all_fastqs/sra
DOWNLOAD_OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/all_fastqs/fastq/$SRR
TMPDIR=/lscratch/$SLURM_JOBID

module load sratoolkit

prefetch $SRR \
    -O $PREFETCH_OUTDIR \
    --max-size 45g \
    --force all
fasterq-dump $PREFETCH_OUTDIR/${SRR}/${SRR}.sra \
    -t $TMPDIR \
    -O $DOWNLOAD_OUTDIR \
    --skip-technical
gzip $DOWNLOAD_OUTDIR/*.fastq
find $DOWNLOAD_OUTDIR -type f -name "*.fastq" -not -name "*.fastq.gz" -delete
rm -r $PREFETCH_OUTDIR/${SRR}
rm -r $TMPDIR/*

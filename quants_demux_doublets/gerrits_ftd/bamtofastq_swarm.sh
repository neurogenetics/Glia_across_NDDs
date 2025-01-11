#! /bin/bash

FASTQ=${1}
BAM_DIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/bams
FASTQ_OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/fastqs

module load cellranger/8.0.1

bamtofastq --nthreads=10 $BAM_DIR/${FASTQ}.bam.1 $FASTQ_OUTDIR/$FASTQ

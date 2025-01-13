#! /bin/bash

FASTQ=${1}

OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/combined
DOUBLETDETECTION_OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/doubletdetection
SCDBLFINDER_OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/scdblfinder
SCDS_OUTDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/scds

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/.singularity

singularity exec --bind /data/ADRD/gerrits_ftd_snRNA/fastq_processing Demuxafy.sif Combine_Results.R \
  -o $OUTDIR/combined_results.tsv \
  -t $DOUBLETDETECTION_OUTDIR \
  -n $SCDBLFINDER_OUTDIR \
  -c $SCDS_OUTDIR 

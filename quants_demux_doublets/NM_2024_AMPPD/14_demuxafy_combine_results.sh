#! /bin/bash

FASTQ=${1}

OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/combined
VIREO_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/vireo
DEMUXALOT_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/demuxalot
DOUBLETDETECTION_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/doubletdetection
SCDBLFINDER_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/scdblfinder
SCDS_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/scds

module load singularity

singularity exec --bind /data/ADRD/amp_pd/transcriptomics/fastq_processing Demuxafy.sif Combine_Results.R \
  -o $OUTDIR/combined_results.tsv \
  -z $DEMUXALOT_OUTDIR \
  -v $VIREO_OUTDIR \
  -t $DOUBLETDETECTION_OUTDIR \
  -n $SCDBLFINDER_OUTDIR \
  -c $SCDS_OUTDIR 

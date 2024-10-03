#!/bin/bash

# need to re-run demuxalot for 3 samples which exceeded memory capacity, only re-running actual demuxalot

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/.singularity

DEMUXALOT_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/demuxalot
BARCODES=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix/barcodes.tsv
BAM=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/possorted_genome_bam.bam
INDS=/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/$FASTQ/donor_list.txt
FILTERED_VCF=/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/$FASTQ/filtered_pool_SNPs.vcf
THREADS=-1 #this along w/ -p in the exec forces it to use the max #CPUs available


singularity exec --bind /data/ADRD/amp_pd/transcriptomics/fastq_processing Demuxafy.sif Demuxalot.py \
    -o $DEMUXALOT_OUTDIR \
    -b $BARCODES \
    -a $BAM \
    -n $INDS \
    -v $FILTERED_VCF \
    -p $THREADS \
    ${CELL_TAG:+-c $CELL_TAG} \
    ${UMI_TAG:+-u $UMI_TAG} \
    -r true

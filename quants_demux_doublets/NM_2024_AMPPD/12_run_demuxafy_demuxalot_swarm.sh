#!/bin/bash

FASTQ=${1}

# filter VCF for only individuals in the pool
INDS=/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/$FASTQ/donor_list.txt
VCF=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/general_inputs/snRNA_SNPs_reordered_tagged_chr.vcf.bgz
FILTERED_VCF=/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/$FASTQ/filtered_pool_SNPs.vcf
BAM=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/possorted_genome_bam.bam

module load bcftools
module load samtools

bcftools view -S $INDS -o $FILTERED_VCF -Ov $VCF

samtools index $BAM

##########

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/.singularity

BARCODES=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix/barcodes.tsv
THREADS=-1 #this along w/ -p in the exec forces it to use the max #CPUs available
DEMUXALOT_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/demuxalot

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

#! /bin/bash

FASTQ=${1}
BARCODES=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix/barcodes.tsv
BAM=/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/possorted_genome_bam.bam
VIREO_OUTDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/vireo
FILTERED_VCF=/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/$FASTQ/filtered_pool_SNPs.vcf
N=6
CELL_TAG=CB
UMI_TAG=UB
FIELD='GT'

module load singularity

export SINGULARITY_CACHEDIR=/data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/.singularity

singularity exec --bind /data/ADRD/amp_pd/transcriptomics/fastq_processing Demuxafy.sif cellsnp_pileup.py \
    -s $BAM \
    -b $BARCODES \
    -O $VIREO_OUTDIR \
    -R $FILTERED_VCF \
    -p 20 \
    --minMAF 0.1 \
    --minCOUNT 20 \
    --cellTAG $CELL_TAG \
    --UMItag $UMI_TAG \
    --gzip true
singularity exec --bind /data/ADRD/amp_pd/transcriptomics/fastq_processing Demuxafy.sif bgzip -c $FILTERED_VCF > $FILTERED_VCF.gz
singularity exec --bind /data/ADRD/amp_pd/transcriptomics/fastq_processing Demuxafy.sif tabix -p vcf $FILTERED_VCF.gz
singularity exec --bind /data/ADRD/amp_pd/transcriptomics/fastq_processing Demuxafy.sif bcftools view $FILTERED_VCF.gz \
    -R $VIREO_OUTDIR/cellSNP.base.vcf.gz \
    -Ov \
    -o $VIREO_OUTDIR/donor_subset.vcf
singularity exec --bind /data/ADRD/amp_pd/transcriptomics/fastq_processing Demuxafy.sif vireo \
    -c $VIREO_OUTDIR \
    -d $VIREO_OUTDIR/donor_subset.vcf \
    -o $VIREO_OUTDIR \
    -t $FIELD

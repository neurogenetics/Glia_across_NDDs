#! /bin/bash
#SBATCH --job-name=nfcore
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --gres=lscratch:800
#SBATCH --time=24:00:00

cd /data/ADRD/glia_across_NDDs

rm -r ./analysis/cellculture/iMG_bulk_RNA/work/*

module load nextflow
export NXF_SINGULARITY_CACHEDIR=/data/ADRD/glia_across_NDDs/nxf_singularity_cache;
export SINGULARITY_CACHEDIR=/data/ADRD/glia_across_NDDs/.singularity;
export TMPDIR=/lscratch/$SLURM_JOB_ID
export NXF_JVM_ARGS="-Xms2g -Xmx4g"

nextflow run nf-core/rnaseq -r 3.19.0 \
-profile biowulf \
-work-dir ./analysis/cellculture/iMG_bulk_RNA/work \
--input ./code_organized/cellculture/iMG_bulk_RNA/fastq_processing/4_batch3_4_DA10682_DA10725/samplesheet_4.csv \
--outdir ./analysis/cellculture/iMG_bulk_RNA/tmp_outs \
--gtf /data/ADRD/human_brain_atlasing/3_fastq_processing/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
--fasta /data/ADRD/human_brain_atlasing/3_fastq_processing/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
--star_index ./analysis/cellculture/iMG_bulk_RNA/reference/index/star \
--salmon_index ./analysis/cellculture/iMG_bulk_RNA/reference/index/salmon \
--trimmer trimgalore \
--aligner star_salmon \
--skip_deseq2_qc true

mv ./analysis/cellculture/iMG_bulk_RNA/tmp_outs/* ./analysis/cellculture/iMG_bulk_RNA/nfcore_outs/batch3_4_DA10682_DA10725/processing_batch4/

#!/bin/bash
#SBATCH --job-name=DESeq2_replication
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --gres=lscratch:5

module load R

Rscript /data/ADRD/glia_across_NDDs/code_organized/microglia/DAM_signature_testing/replication_organized/3.1_run_DESeq2_replication_data.R
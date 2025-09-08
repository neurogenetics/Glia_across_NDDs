#!/bin/bash
#SBATCH --job-name=dream_analysis
#SBATCH --output=logs/dream_%A_%a.out
#SBATCH --error=logs/dream_%A_%a.err
#SBATCH --mem=64g                 
#SBATCH --cpus-per-task=24         
#SBATCH --gres=lscratch:5
#SBATCH --time=04:00:00            
#SBATCH --array=0-3

mkdir -p logs

batches=("XR10585" "DA10634" "DA10682" "DA10725")
current_batch=${batches[$SLURM_ARRAY_TASK_ID]}

echo "Starting job $SLURM_JOB_ID, task $SLURM_ARRAY_TASK_ID, processing batch: $current_batch"

module load R

# Run your modified R script, passing the current batch ID as an argument
Rscript job_run_dream_batch.R $current_batch

echo "Finished task $SLURM_ARRAY_TASK_ID."
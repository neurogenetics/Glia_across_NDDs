#!/bin/bash
#SBATCH --job-name=cNMF_microglia
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=20G

# this was run as a batch job!

export PATH=$HOME/.local/bin:$PATH

DIR=/data/ADRD/glia_across_NDDs/analysis/microglia/cluster_characterization/cNMF_by_dataset/Pineda

~/.local/bin/cnmf prepare --output-dir $DIR --name microglia -c $DIR/microglia.Corrected.HVG.Varnorm.h5ad --tpm $DIR/microglia.TP10K.h5ad --genes-file $DIR/microglia.Corrected.HVGs.txt -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 -n 100 --seed 14 --numgenes 2000 

~/.local/bin/cnmf factorize --output-dir $DIR --name microglia --worker-index 0 --total-workers 1

~/.local/bin/cnmf combine --output-dir $DIR --name microglia

~/.local/bin/cnmf k_selection_plot --output-dir $DIR --name microglia

########################################
########################################
########################################

# this was run in interactive session

DIR=/data/ADRD/glia_across_NDDs/analysis/microglia/cluster_characterization/cNMF_by_dataset/Pineda

# run once w/ high density threshold
~/.local/bin/cnmf consensus --output-dir $DIR --name microglia --components 12 --local-density-threshold 2 --show-clustering

# then run again after checking density threshold plot
~/.local/bin/cnmf consensus --output-dir $DIR --name microglia --components 12 --local-density-threshold 0.10 --show-clustering

# bamtofastq (cellranger)
swarm -f /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts/bamtofastq_swarm.swarm \
    -g 50 \
    -t 10 \
    --module cellranger \
    --logdir /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs \
    --time=24:00:00 \
    --job-name swarm_bamtofastq_gerrits

# cellranger
swarm -f /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts/cellranger_swarm.swarm \
    -g 124 \
    -t 24 \
    --module cellranger \
    --logdir /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs \
    --time=18:00:00 \
    --gres=lscratch:700 \
    --job-name swarm_cellranger_gerrits

# cellbender
swarm -f /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts/cellbender_swarm.swarm \
    -g 100 \
    -t 12 \
    --module cellbender \
    --logdir /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs \
    --time=12:00:00 \
    --partition=gpu \
    --gres=gpu:k80:1 \
    --job-name swarm_cellbender_gerrits

# doubletdetection
swarm -f /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts/doubletdetection_swarm.swarm \
    -g 50 \
    -t 12 \
    -b 2 \
    --module singularity \
    --logdir /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs \
    --time=4:00:00 \
    --job-name swarm_doubletdetection_gerrits

# scdblfinder
swarm -f /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts/scdblfinder_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs \
    --time=2:00:00 \
    --job-name swarm_scdblfinder_gerrits

# scds
swarm -f /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts/scds_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs \
    --time=2:00:00 \
    --job-name swarm_scds_gerrits

# combine results
swarm -f /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts/demuxafy_combine_results.swarm \
    -g 3 \
    -t 1 \
    -b 20 \
    --module singularity \
    --logdir /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs \
    --time=00:10:00 \
    --job-name swarm_combine_results_gerrits

# filter doublets
swarm -f /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts/filter_doublets.swarm \
    -g 50 \
    -t 1 \
    --module R \
    --logdir /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs \
    --time=2:00:00 \
    --job-name swarm_filter_cells_gerrits \
    --gres=lscratch:5

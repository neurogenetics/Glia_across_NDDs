# cellranger
swarm -f /data/ADRD/AD_multiregion/fastq_processing/swarmscripts/cellranger_swarm.swarm \
    -g 124 \
    -t 24 \
    --module cellranger \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=18:00:00 \
    --gres=lscratch:700 \
    --job-name swarm_cellranger_ad

# cellbender
swarm -f /data/ADRD/AD_multiregion/fastq_processing/swarmscripts/cellbender_swarm.swarm \
    -g 100 \
    -t 12 \
    -b 4 \
    --module cellbender \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=12:00:00 \
    --partition=gpu \
    --gres=gpu:k80:1 \
    --job-name swarm_cellbender_ad

# doubletdetection
swarm -f /data/ADRD/AD_multiregion/fastq_processing/swarmscripts/doubletdetection_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=4:00:00 \
    -b 5 \
    --job-name swarm_doubletdetection_ad

# scdblfinder
swarm -f /data/ADRD/AD_multiregion/fastq_processing/swarmscripts/scdblfinder_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=2:00:00 \
    -b 10 \
    --job-name swarm_scdblfinder_ad

# scds
swarm -f /data/ADRD/AD_multiregion/fastq_processing/swarmscripts/scds_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=2:00:00 \
    -b 10 \
    --job-name swarm_scds_ad

# combine results
swarm -f /data/ADRD/AD_multiregion/fastq_processing/swarmscripts/demuxafy_combine_results.swarm \
    -g 3 \
    -t 1 \
    -b 20 \
    --module singularity \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=00:10:00 \
    --job-name swarm_combine_results_ad

# filter doublets
swarm -f /data/ADRD/AD_multiregion/fastq_processing/swarmscripts/filter_doublets.swarm \
    -g 50 \
    -t 1 \
    -b 10 \
    --module R \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=2:00:00 \
    --job-name swarm_filter_cells_ad \
    --gres=lscratch:5

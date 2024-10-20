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
    -b 2 \
    --module cellbender \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=12:00:00 \
    --partition=gpu \
    --gres=gpu:k80:1 \
    --job-name swarm_cellbender_alsftd

# sratoolkit (prefetch + fasterq-dump)
swarm -f /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmscripts/prefetch_fasterq_dump_swarm.swarm \
    -g 20 \
    -t 10 \
    -b 2 \
    --module sratoolkit \
    --logdir /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmlogs \
    --time=12:00:00 \
    --gres=lscratch:500 \
    --job-name swarm_sra_alsftd

# cellranger
swarm -f /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmscripts/cellranger_swarm.swarm \
    -g 124 \
    -t 24 \
    --module cellranger \
    --logdir /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmlogs \
    --time=18:00:00 \
    --gres=lscratch:700 \
    --job-name swarm_cellranger_alsftd

# cellbender
swarm -f /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmscripts/cellbender_swarm.swarm \
    -g 100 \
    -t 12 \
    -b 2 \
    --module cellbender \
    --logdir /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmlogs \
    --time=12:00:00 \
    --partition=gpu \
    --gres=gpu:k80:1 \
    --job-name swarm_cellbender_alsftd

# doubletdetection
swarm -f /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmscripts/doubletdetection_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmlogs \
    --time=4:00:00 \
    -b 3 \
    --job-name swarm_doubletdetection_alsftd

# scdblfinder
swarm -f /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmscripts/scdblfinder_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmlogs \
    --time=2:00:00 \
    -b 10 \
    --job-name swarm_scdblfinder_alsftd

# scds
swarm -f /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmscripts/scds_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmlogs \
    --time=2:00:00 \
    -b 5\
    --job-name swarm_scds_alsftd

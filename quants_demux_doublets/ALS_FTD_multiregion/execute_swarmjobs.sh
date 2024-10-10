# sratoolkit (prefetch + fasterq-dump)
swarm -f /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmscripts/prefetch_fasterq_dump_swarm.swarm \
    -g 10 \
    -t 10 \
    -b 20 \
    --module sratoolkit \
    --logdir /data/ADRD/ALSFTD_multiregion/fastq_processing/swarmlogs \
    --time=4:00:00 \
    --gres=lscratch:500 \
    --job-name swarm_sra_alsftd

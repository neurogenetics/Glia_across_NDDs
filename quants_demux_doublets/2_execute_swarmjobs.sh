# commands for executing swarm scripts made in 1_make_swarmscripts.sh

# cellranger
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/cellranger_swarm.swarm \
    -g 124 \
    -t 24 \
    --module cellranger \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=36:00:00 \
    --gres=lscratch:700 \
    --job-name swarm_cellranger_amppd
# cellbender
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/cellbender_swarm.swarm \
    -g 100 \
    -t 12 \
    --module cellbender \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=24:00:00 \
    --partition=gpu \
    --gres=gpu:k80:1 \
    --job-name swarm_cellbender_amppd
# doubletdetection
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/doubletdetection_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=24:00:00 \
    --job-name swarm_doubletdetection_amppd
# scdblfinder
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/scdblfinder_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=24:00:00 \
    --job-name swarm_scdblfinder_amppd
# scds
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/scds_swarm.swarm \
    -g 50 \
    -t 12 \
    --module singularity \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=24:00:00 \
    --job-name swarm_scds_amppd
# demuxalot
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/demuxalot_swarm.swarm \
    -g 100 \
    -t 24 \
    --module singularity,bcftools,samtools \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=24:00:00 \
    --job-name swarm_demuxalot_amppd
# vireo
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/vireo_swarm.swarm \
    -g 20 \
    -t 24 \
    --module singularity \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=48:00:00 \
    --job-name swarm_vireo_amppd
# combine results
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/demuxafy_combine_results.swarm \
    -g 3 \
    -t 1 \
    --module singularity \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=00:10:00 \
    --job-name swarm_combine_results_amppd
# filter doublets
swarm -f /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/filter_doublets.swarm \
    -g 50 \
    -t 1 \
    --module R \
    --logdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmlogs \
    --time=2:00:00 \
    --job-name swarm_filter_cells_amppd \
    --gres=lscratch:5

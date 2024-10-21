# had to re-run 3 files bc gpu node failed

SWARM_FILENAME=/data/ADRD/AD_multiregion/fastq_processing/swarmscripts/cellbender_rerun_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/AD_multiregion/fastq_processing/cellbender_rerun.txt)
JOB_SCRIPT=/data/ADRD/AD_multiregion/fastq_processing/jobscripts/run_cellbender_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

swarm -f /data/ADRD/AD_multiregion/fastq_processing/swarmscripts/cellbender_rerun_swarm.swarm \
    -g 100 \
    -t 12 \
    --module cellbender \
    --logdir /data/ADRD/AD_multiregion/fastq_processing/swarmlogs \
    --time=12:00:00 \
    --partition=gpu \
    --gres=gpu:k80:1 \
    --job-name swarm_cellbender_rerun

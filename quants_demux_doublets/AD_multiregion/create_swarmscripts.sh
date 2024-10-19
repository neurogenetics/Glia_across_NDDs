SWARM_FILENAME=/data/ADRD/AD_multiregion/fastq_processing/swarmscripts/cellranger_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/AD_multiregion/fastq_processing/split_fastqs.txt)
JOB_SCRIPT=/data/ADRD/AD_multiregion/fastq_processing/jobscripts/run_cellranger_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

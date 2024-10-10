SWARM_FILENAME=/data/ADRD/ALSFTD_multiregion/fastq_processing/swarmscripts/prefetch_fasterq_dump_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/ALSFTD_multiregion/fastq_processing/cleaned_SRRs.txt)
JOB_SCRIPT=/data/ADRD/ALSFTD_multiregion/fastq_processing/jobscripts/prefetch_fasterq_dump_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

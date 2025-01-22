SWARM_FILENAME=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/swarmscripts/prefetch_fasterq_dump_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/Drager_SRR.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/jobscripts/prefetch_fasterq_dump_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/swarmscripts/cellranger_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/Drager_SRR.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/jobscripts/run_cellranger_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/swarmscripts/cellbender_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/Drager_SRR.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/jobscripts/run_cellbender_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/swarmscripts/doubletdetection_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/Drager_SRR.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/jobscripts/run_demuxafy_doubletdetection_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/swarmscripts/scdblfinder_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/Drager_SRR.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/jobscripts/run_demuxafy_scdblfinder_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/swarmscripts/scds_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/Drager_SRR.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/jobscripts/run_demuxafy_scds_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/swarmscripts/demuxafy_combine_results.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/Drager_SRR.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/jobscripts/demuxafy_combine_results.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/swarmscripts/filter_doublets.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/Drager_SRR.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/fastq_processing/jobscripts/filter_doublets.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done


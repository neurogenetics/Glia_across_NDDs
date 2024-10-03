# bash script to create swarm files, run in command-line

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/cellranger_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/4_run_cellranger_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/cellbender_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/fastq_list_rerun_demuxalot.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/8_run_cellbender_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/doubletdetection_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/9_run_demuxafy_doubletdetection_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/scdblfinder_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/10_run_demuxafy_scdblfinder_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/scds_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/11_run_demuxafy_scds_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/demuxalot_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/12_run_demuxafy_demuxalot_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/demuxalot_rerun_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/demuxalot_rerun_highmem.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/vireo_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/13_run_demuxafy_vireo_swarm.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/demuxafy_combine_results.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/14_demuxafy_combine_results.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

##########

SWARM_FILENAME=/data/ADRD/amp_pd/transcriptomics/fastq_processing/swarmscripts/filter_doublets.swarm
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)
JOB_SCRIPT=/data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/15_filter_doublets.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done

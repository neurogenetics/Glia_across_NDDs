# create swarmscript

SWARM_FILENAME=/data/ADRD/glia_across_NDDs/code_organized/microglia/differential_expression/1.4_run_microglia_diffexp_swarm.swarm
FASTQ_LIST=$(cat /data/ADRD/glia_across_NDDs/analysis/diffexp_tests.txt)
JOB_SCRIPT=/data/ADRD/glia_across_NDDs/code_organized/microglia/differential_expression/1.2_run_microglia_diffexp_bash.sh

for FASTQ in ${FASTQ_LIST}
do
echo ${JOB_SCRIPT} ${FASTQ} >> ${SWARM_FILENAME}
done


# execute swarm job

swarm -f /data/ADRD/glia_across_NDDs/code_organized/microglia/differential_expression/1.4_run_microglia_diffexp_swarm.swarm \
    -g 100 \
    -t 2 \
    --module R \
    --logdir /data/ADRD/glia_across_NDDs/analysis/swarmlogs \
    --time=24:00:00 \
    --gres=lscratch:5 \
    --job-name microglia_diffexp
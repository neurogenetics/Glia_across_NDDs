# Drager 2022
For this analysis work_dir = /data/ADRD/glia_across_NDDs/iMG_meta_analysis/Drager2022/

**Workflow for code:**
## 0. Set up environment
- 0A - create a dir "{work_dir}/fastq_processing"
   - 0B - create an "{work_dir}/fastq_processing/Drager_SRR.txt" containing all SRRs deliminaited by line 
 
## 1. Initiliaze Directories and Files
- 1A - Make output directories "{work_dir}/make_output_dirs.sh"
   - This step will use the file in fastq_processing to create necessary directories
```shell
{work_dir}/make_output_dirs.sh
#If there is an error, make sure this file is executable
#optional (one file):
# chmod 777 {work_dir}/make_output_dirs.sh 
#optional (all files in dir):
# chmod -R 777 {work_dir}

```
- 1B - Make swarm scripts "{work_dir}/create_swarmscripts.sh"
   - Change the text so that SWARM_FILNAME = {work_dir}/fastq_processing/swarmscripts/*
   - Ensure FASTQ_LIST = {work_dir}/fastq_processing/*_SRR.txt
   - Change the text so that JOB_SCRIPT = {work_dir}/fastq_processing/jobscripts/*
```shell
{work_dir}/create_swarmscripts.sh

```
 - 1C - Download fastqs from SRA as swarm 
   - file 1, copy from other project - {work_dir}/fastq_processing/jobscripts/prefetch_fasterq_dump_swarm.sh
   - file 2, created by step 1B - {work_dir}/fastq_processing/swarmscripts/prefetch_fasterq_dump_swarm.swarm
`
```shell
swarm -f {work_dir}/fastq_processing/swarmscripts/prefetch_fasterq_dump_swarm.swarm \
    -g 20 \
    -t 10 \
    -b 2 \
    --module sratoolkit \
    --logdir {work_dir}/fastq_processing/swarmlogs \
    --time=12:00:00 \
    --gres=lscratch:500 \
    --job-name swarm_sra_drager


```

 - 1D - Rename fastqs --> differs based on fastq names, needs to split, or rename
   - eg renaming: https://github.com/neurogenetics/DLB_VHD/blob/main/building_EC_reference/fastq_processing/franjic_2022/rename_fastqs.sh
   - however Drager R1/R2 are separate SRR's so I decided to do manually: see {work_dir}/Drager_SraRunTable(1).csv/{'Run'|'Sample Name'}
   - Selected naming for this process: {Sample Name}_S1_{R1 = lowerSRR, R2 = higher SRR}_001.fastq.gz
 
```shell
mkdir {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387652
mv {work_dir}/fastq_processing/all_fastqs/fastq/SRR14828083/SRR14828083.fastq.gz {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387652/GSM5387652_S1_R1_001.fastq.gz
mv {work_dir}/fastq_processing/all_fastqs/fastq/SRR14828084/SRR14828084.fastq.gz {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387652/GSM5387652_S1_R2_001.fastq.gz

mkdir {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387653
mv {work_dir}/fastq_processing/all_fastqs/fastq/SRR14828085/SRR14828085.fastq.gz {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387653/GSM5387653_S1_R1_001.fastq.gz
mv {work_dir}/fastq_processing/all_fastqs/fastq/SRR14828086/SRR14828086.fastq.gz {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387653/GSM5387653_S1_R2_001.fastq.gz

mkdir {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387654
mv {work_dir}/fastq_processing/all_fastqs/fastq/SRR14828087/SRR14828087.fastq.gz {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387654/GSM5387654_S1_R1_001.fastq.gz
mv {work_dir}/fastq_processing/all_fastqs/fastq/SRR14828088/SRR14828088.fastq.gz {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387654/GSM5387654_S1_R2_001.fastq.gz

mkdir {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387655
mv {work_dir}/fastq_processing/all_fastqs/fastq/SRR14828089/SRR14828089.fastq.gz {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387655/GSM5387655_S1_R1_001.fastq.gz
mv {work_dir}/fastq_processing/all_fastqs/fastq/SRR14828090/SRR14828090.fastq.gz {work_dir}/fastq_processing/all_fastqs/fastq/GSM5387655/GSM5387655_S1_R2_001.fastq.gz
```
_now in the swarm files pointing to SRR numbers need to change them to GSM numbers, will make sure to note where that happened._


## 2. Cellranger
 - alter {work_dir}/swarmscripts/cellranger_swarm.swarm to GSM instead of SRR, this is an exception because of how data was inputted. no alterations to {work_dir}/jobscripts/run_cellranger_swarm.sh needed
```bash
swarm -f {work_dir}/fastq_processing/swarmscripts/cellranger_swarm.swarm \
    -g 124 \
    -t 24 \
    --module cellranger \
    --logdir {work_dir}/fastq_processing/swarmlogs \
    --time=18:00:00 \
    --gres=lscratch:700 \
    --job-name swarm_cellranger_drager

```

## 3. Ambient RNA

## 4. Ambient RNA

## 5. Filtering in Seurat


2. Cell-level QC, generating microglia/astrocyte subsets for each dataset (**qc_filtering**)

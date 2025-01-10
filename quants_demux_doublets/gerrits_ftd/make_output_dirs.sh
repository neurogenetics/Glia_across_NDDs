dos2unix gerrits_srrs.txt

FASTQ_LIST=$(cat /data/ADRD/gerrits_ftd_snRNA/fastq_processing/gerrits_srrs.txt)

for FASTQ in ${FASTQ_LIST}
do
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmscripts
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/bams
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/swarmlogs
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/jobscripts
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/all_fastqs/fastq/$FASTQ
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/all_fastqs/sra
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/cellranger/outs/$FASTQ
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/cellbender/outs/$FASTQ
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/doubletdetection
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/scdblfinder
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/scds
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/demuxafy/outs/$FASTQ/combined
mkdir -p /data/ADRD/gerrits_ftd_snRNA/fastq_processing/final_outs/$FASTQ
done

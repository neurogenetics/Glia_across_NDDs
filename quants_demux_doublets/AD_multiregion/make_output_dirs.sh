FASTQ_LIST=$(cat /data/ADRD/AD_multiregion/fastq_processing/split_fastqs.txt)

for FASTQ in ${FASTQ_LIST}
do
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/cellranger/outs/$FASTQ
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/cellbender/outs/$FASTQ
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/doubletdetection
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/scdblfinder
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/scds
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/combined
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/final_outs/$FASTQ
done

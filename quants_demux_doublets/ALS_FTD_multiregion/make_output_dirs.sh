FASTQ_LIST=$(cat /data/ADRD/ALSFTD_multiregion/fastq_processing/cleaned_SRRs.txt)

for FASTQ in ${FASTQ_LIST}
do
mkdir /data/ADRD/ALSFTD_multiregion/fastq_processing/all_fastqs/fastq/$FASTQ
mkdir /data/ADRD/ALSFTD_multiregion/fastq_processing/cellranger/outs/$FASTQ
mkdir /data/ADRD/ALSFTD_multiregion/fastq_processing/cellbender/outs/$FASTQ
mkdir /data/ADRD/ALSFTD_multiregion/fastq_processing/demuxafy/outs/$FASTQ
mkdir /data/ADRD/ALSFTD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/doubletdetection
mkdir /data/ADRD/ALSFTD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/scdblfinder
mkdir /data/ADRD/ALSFTD_multiregion/fastq_processing/demuxafy/outs/$FASTQ/scds
done

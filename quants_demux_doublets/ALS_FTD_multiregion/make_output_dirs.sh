FASTQ_LIST=$(cat /data/ADRD/ALSFTD_multiregion/fastq_processing/all_fastqs/ALS_SRRs.txt)

for FASTQ in ${FASTQ_LIST}
do
mkdir /data/ADRD/ALSFTD_multiregion/fastq_processing/all_fastqs/fastq/$FASTQ
done

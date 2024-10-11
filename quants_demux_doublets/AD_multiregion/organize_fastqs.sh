FASTQ_LIST=$(cat /data/ADRD/AD_multiregion/fastq_processing/cleaned_AD_fastqs.txt)

for FASTQ in ${FASTQ_LIST}
do
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/all_fastqs/$FASTQ
mv /data/ADRD/AD_multiregion/fastqs/'Sequencing Round 1'/*$FASTQ /data/ADRD/AD_multiregion/fastq_processing/all_fastqs/$FASTQ/
mv /data/ADRD/AD_multiregion/fastqs/'Sequencing Round 2'/*$FASTQ /data/ADRD/AD_multiregion/fastq_processing/all_fastqs/$FASTQ/
done

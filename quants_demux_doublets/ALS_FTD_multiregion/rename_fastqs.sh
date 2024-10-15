FASTQ_LIST=$(cat /data/ADRD/ALSFTD_multiregion/fastq_processing/cleaned_SRRs.txt)

for FASTQ in FASTQ_LIST
do
FASTQ_DIR=/data/ADRD/ALSFTD_multiregion/fastq_processing/all_fastqs/fastq/$FASTQ
mv ${FASTQ_DIR}/${FASTQ}_1.fastq.gz ${FASTQ_DIR}/${FASTQ}_S1_R1_001.fastq.gz
mv ${FASTQ_DIR}/${FASTQ}_2.fastq.gz ${FASTQ_DIR}/${FASTQ}_S1_R2_001.fastq.gz
done

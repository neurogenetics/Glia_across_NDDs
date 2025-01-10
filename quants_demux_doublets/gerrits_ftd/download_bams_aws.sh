FASTQ_LIST=$(cat /data/ADRD/gerrits_ftd_snRNA/fastq_processing/gerrits_srrs.txt)

for FASTQ in ${FASTQ_LIST}
do
aws --no-sign-request s3 cp s3://sra-pub-src-2/$FASTQ/ ./bams --recursive
echo "$FASTQ done"
done

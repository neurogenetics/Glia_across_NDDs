# run in helix instance
# configure aws w/ region as us-east-1

FASTQ_LIST=$(cat /data/ADRD/gerrits_ftd_snRNA/fastq_processing/gerrits_srrs.txt)

for FASTQ in ${FASTQ_LIST}
do
aws --no-sign-request s3 cp s3://sra-pub-src-2/$FASTQ/ ./bams --recursive
echo "$FASTQ done"
done

# for some reason some were in a different bucket, repeat w/ different bucket name

FASTQ_LIST=$(cat /data/ADRD/gerrits_ftd_snRNA/fastq_processing/gerrits_srrs.txt)

for FASTQ in ${FASTQ_LIST}
do
aws --no-sign-request s3 cp s3://sra-pub-src-1/$FASTQ/ ./bams --recursive
echo "$FASTQ done"
done

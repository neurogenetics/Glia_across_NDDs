# run in command-line, HELIX!!!

FASTQ_LIST=$(cat /data/ADRD/AD_multiregion/fastq_processing/round_1_fastqs.txt)

for FASTQ in ${FASTQ_LIST}
do
gunzip --keep /data/ADRD/AD_multiregion/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix/barcodes.tsv.gz
done

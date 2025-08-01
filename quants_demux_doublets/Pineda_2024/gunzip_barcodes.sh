# run in command-line, HELIX!!!

FASTQ_LIST=$(cat /data/ADRD/ALSFTD_multiregion/fastq_processing/cleaned_SRRs.txt)

for FASTQ in ${FASTQ_LIST}
do
gunzip --keep /data/ADRD/ALSFTD_multiregion/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix/barcodes.tsv.gz
done

# run in command-line, HELIX!!!

FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)

for FASTQ in ${FASTQ_LIST}
do
gunzip --keep /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/filtered_feature_bc_matrix/barcodes.tsv.gz
done

# run in command-line, HELIX!!!

FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/fastq_list.txt)

for FASTQ in ${FASTQ_LIST}
do
gunzip --keep /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/outs/$FASTQ/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
done

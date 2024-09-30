# renaming fastqs to match CellRanger naming conventions
# run in command-line

FASTQ_DIR="/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs"
FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)

for SAMPLE in $FASTQ_LIST; do
    SAMPLE_DIR="$FASTQ_DIR/$SAMPLE"

    for FILE in "$SAMPLE_DIR"/*.fastq.gz; do
        # extract the lane and read type (R1 or R2)
        LANE=$(echo "$FILE" | grep -oP "L\d{3}")
        READ_TYPE=$(echo "$FILE" | grep -oP "R[12]")

        # create the new filename in the CellRanger convention
        NEW_NAME="${SAMPLE}_S1_${LANE}_${READ_TYPE}_001.fastq.gz"

        # rename the file
        mv "$FILE" "$SAMPLE_DIR/$NEW_NAME"
    done
done

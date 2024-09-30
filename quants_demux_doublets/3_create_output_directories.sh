# run in command-line

FASTQ_LIST=$(cat /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt)

for FASTQ in ${FASTQ_LIST}
do
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/doubletdetection
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/scdblfinder
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/scds
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/demuxalot
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/demuxlet
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/vireo
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/demuxafy/samples_outs/$FASTQ/combined
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/cellbender/outs/$FASTQ
mkdir /data/ADRD/amp_pd/transcriptomics/fastq_processing/final_outs/$FASTQ
done

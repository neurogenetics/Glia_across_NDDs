# THIS IS JUST HOW I DID IT FOR THE AMP-PD SAMPLES!!!
# NAMING CONVENTIONS WILL BE DIFFERENT!!!
# Essentially what you need is a .txt file containing the names of the samples, exactly as they appear in the VCF, separated by a line

fastq_list <- read.delim('/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt', header = F)
fastqs <- fastq_list$V1

for (fastq in fastqs){
  txt_directory <- paste0('/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/', fastq, '/', fastq, '.provided_sample_list.txt')
  txt <- read.delim(txt_directory)
  txt$sample_id <- gsub('-BLM0-[A-Z]{3,4}-RSN', '', txt$sample_id)
  new_txt_directory <- paste0('/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/', fastq, '/donor_list.txt')
  write.table(txt$sample_id, file = new_txt_directory, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

for (fastq in fastqs){
  txt_directory <- paste0('/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/', fastq, '/', fastq, '.provided_sample_list.txt')
  txt <- read.delim(txt_directory)
  txt$sample_id <- gsub('-BLM0', '', txt$sample_id)
  txt$sample_id <- gsub('-RSN', '', txt$sample_id)
  new_txt_directory <- paste0('/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/', fastq, '/donor_region_list.txt')
  write.table(txt$sample_id, file = new_txt_directory, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

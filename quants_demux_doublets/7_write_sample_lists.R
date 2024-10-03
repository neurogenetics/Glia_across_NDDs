fastq_list <- read.delim('/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt', header = F)
fastqs <- fastq_list$V1

for (fastq in fastqs){
  txt_directory <- paste0('/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/', fastq, '/', fastq, '.provided_sample_list.txt')
  txt <- read.delim(txt_directory)
  txt$sample_id <- gsub('-BLM0-[A-Z]{3,4}-RSN', '', txt$sample_id)
  new_txt_directory <- paste0('/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/', fastq, '/donor_list.txt')
  write.table(txt$sample_id, file = new_txt_directory, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# !!! needed to manually edit 3 sets (for both duplicates): Set4, Set5, Set7a (both duplicates, so 6 total pools)
# they loaded 2 brain regions from the same donor in these pools, which can't be separated
# cells from the dupe donor will be discarded
# just edited the sample_id txt file to remove one of the copies of that donor

for (fastq in fastqs){
  txt_directory <- paste0('/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/', fastq, '/', fastq, '.provided_sample_list.txt')
  txt <- read.delim(txt_directory)
  txt$sample_id <- gsub('-BLM0', '', txt$sample_id)
  txt$sample_id <- gsub('-RSN', '', txt$sample_id)
  new_txt_directory <- paste0('/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs/', fastq, '/donor_region_list.txt')
  write.table(txt$sample_id, file = new_txt_directory, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# for those same pools, also manually removed one of the regions from the region txt
# will manually remove all cells from those donors from those pools after combining everything

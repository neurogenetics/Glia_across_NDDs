fastq_directory <- "/data/ADRD/amp_pd/transcriptomics/fastq_processing/all_fastqs"

fastq_names <- list.dirs(fastq_directory, full.names = FALSE, recursive = FALSE)

writeLines(fastq_names, con = "/data/ADRD/amp_pd/transcriptomics/fastq_processing/cellranger/all_fastqs.txt")

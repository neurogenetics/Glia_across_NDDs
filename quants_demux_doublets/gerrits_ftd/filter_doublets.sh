#!/bin/bash

FASTQ=${1}

Rscript /data/ADRD/amp_pd/transcriptomics/fastq_processing/jobscripts/filter_doublets.R $FASTQ

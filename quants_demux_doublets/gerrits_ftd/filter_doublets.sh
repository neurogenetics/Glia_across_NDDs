#!/bin/bash

FASTQ=${1}

Rscript /data/ADRD/gerrits_ftd_snRNA/fastq_processing/jobscripts/filter_doublets.R $FASTQ

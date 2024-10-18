#!/bin/bash

FASTQ=${1}

Rscript /data/ADRD/ALSFTD_multiregion/fastq_processing/jobscripts/filter_doublets.R $FASTQ

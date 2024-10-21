#!/bin/bash

FASTQ=${1}

Rscript /data/ADRD/AD_multiregion/fastq_processing/jobscripts/filter_doublets.R $FASTQ

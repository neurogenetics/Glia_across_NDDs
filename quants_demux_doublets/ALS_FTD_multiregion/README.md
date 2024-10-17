Samples are not multiplexed! Only need to do general doublet detection
* Create and execute swarm scripts (**create_swarmscripts.sh and execute_swarmjobs.sh**)
1. Make output directories (**make_output_dirs.sh**)
2. Download fastqs from SRA (**prefetch_fasterq_dump.sh**)
3. Rename fastqs to match CellRanger naming convention (**rename_fastqs.sh**)
4. Run CellRanger (**run_cellranger_swarm.sh**)
5. Unzip barcodes.tsv.gz for Demuxafy (**gunzip_barcodes.sh**)
6. Run CellBender (**run_cellbender_swarm.sh**)

*Samples are not multiplexed! Only need to do general doublet detection*

* Create and execute swarm scripts (create_swarmscripts.sh and execute_swarmjobs.sh)

1. Organize fastqs into folders (**organize_fastqs.sh**)
2. Make output directories (**make_output_dirs.sh**)
3. Run CellRanger (**run_cellranger_swarm.sh**)
4. Unzip barcodes.tsv.gz for Demuxafy (**gunzip_barcodes.sh**)
5. Run CellBender (**run_cellbender_swarm.sh**)
6. Run Demuxafy doublet detection packages (**run_demuxafy_doubletdetection_swarm.sh**, **run_demuxafy_scdblfinder_swarm.sh**, and **run_demuxafy_scds_swarm.sh**)
7. Combine results from Demuxafy packages (**demuxafy_combine_results.sh**)
8. Filter doublets and make/save Seurat objects for each sample (**filter_doublets.R** and **filter_doublets.sh**)


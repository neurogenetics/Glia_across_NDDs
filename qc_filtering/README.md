Cells were filtered using the following metrics, calculated from the **raw CellRanger counts** for each dataset separately: 
* 2.5th percentile < # transcripts > 97.5th percentile
* 2.5th percentile < # unique genes
* 2% > % mitochondrial genes
* 2% > % ribosomal genes

1. Get QC metrics from all 3 unfiltered datasets (**get_qc_metrics.R**)
2. Filter data based on UMIs and features (**all_data_qc_filters.R**)
3. Apply QC filters to Seurat objects (**apply_qc_filters.R**)
4. Split each dataset by brain region (**split_by_region.R**)

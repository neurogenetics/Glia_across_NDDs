# Glia across NDDs
*Analysis of three atlas-scale snRNA-seq datasets from multiple neurodegenerative diseases*

Parkinson's disease (AMP-PD, unpublished)

ALS/FTLD (Pineda et al., 2024)

Alzheimer's disease (Mathys et al., 2024)

**Workflow for code:**
1. Download data, generate counts matrices, doublet detection, demultiplexing (**quants_demux_doublets**)
2. Cell-level QC, generating microglia/astrocyte subsets for each dataset (**qc_filtering**)

# Glia across neurodegenerative diseases
*Analysis of four atlas-scale snRNA-seq datasets from multiple neurodegenerative diseases*
## Discovery datasets
ALS/FTLD (Pineda et al., 2024)
Alzheimer's disease (Mathys et al., 2024)
FTD-GRN (Gerrits et al.,)
Parkinson's disease (NM et al., 2024)

## Other analyses
Replication datasets - 11 sample series, pulled as processed data
iPSC-modeling - iMicroglia data produced from 5 PPMI lines (available on PPMI database)

# Explanation of Code
**Workflow for code:**
1. Download data, generate counts matrices, doublet detection, demultiplexing (**quants_demux_doublets**)
2. Cell-level QC, generating microglia/astrocyte subsets for each dataset (**qc_filtering**)
3. Proportions, differential expression, and latent factorization (**analysis**)

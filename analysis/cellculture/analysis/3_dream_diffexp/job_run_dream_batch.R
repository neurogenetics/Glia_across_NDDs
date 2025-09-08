#!/usr/bin/env Rscript
# run_dream_batch.R


# --- Setup: Load libraries and data ---
library(variancePartition)
library(tidyverse)
library(edgeR)
library(BiocParallel) # Explicitly load for bplapply

setwd("/data/ADRD/glia_across_NDDs")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No batch ID provided. Usage: Rscript run_dream_batch.R <batch_id>", call. = FALSE)
}
batch <- args[1]

cts <- readRDS("./analysis/cellculture/iMG_bulk_RNA/cleaned_counts/all_batches_cts_merged.rds")
meta <- readRDS("./analysis/cellculture/iMG_bulk_RNA/metadata/sample_meta_merged.rds")

# --- Parallel processing setup ---
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 4))
print(paste("Setting up parallel backend with", n_cores, "cores."))
param <- SnowParam(workers = n_cores, type = "SOCK", progressbar = TRUE)
register(param)

# --- Analysis: Process the single batch provided ---
print(paste0("Beginning batch: ", batch))
meta_filtered <- meta[grepl(batch, rownames(meta)), ]
cts_filtered <- cts[, colnames(cts) %in% rownames(meta_filtered)]

dge <- DGEList(cts_filtered)
dge <- calcNormFactors(dge, method = "TMM")
cpm <- as.data.frame(edgeR::cpm(dge))

meta_filtered <- meta_filtered[match(colnames(cpm), rownames(meta_filtered)), ]
if("CTL_0.00E+00"%in%meta_filtered$treatment_dose){
  meta_filtered <- meta_filtered%>%
    mutate(treatment_dose = case_when(treatment_dose=="CTL_0.00E+00"~"CTL_0",
                                      treatment_dose!="CTL_0.00E+00"~treatment_dose))
}

genes_keep <- vector("logical", nrow(cpm))
for (t in unique(meta_filtered$treatment)) {
  samples_in_group <- rownames(meta_filtered)[meta_filtered$treatment == t]
  cpm_sub <- cpm[, samples_in_group, drop = FALSE]
  gene_pass <- rowSums(cpm_sub > 1) >= 3
  genes_keep <- genes_keep | gene_pass
}
cts_filtered <- cts_filtered[genes_keep, ]

treatments <- unique(meta_filtered$treatment_dose)
treatments <- treatments[!treatments %in% "CTL_0"]

# --- NEW: Define a function to process ONE treatment ---
# MODIFIED LINE 1: Add the data frames as arguments to the function
run_dream_for_treatment <- function(treatment, local_meta, local_cts) {
  print(paste0("Processing control vs. ", treatment))
  
  meta_test <- local_meta[local_meta$treatment_dose %in% c("CTL_0", treatment), ]
  cts_test <- local_cts[, colnames(local_cts) %in% rownames(meta_test)]
  meta_test$PPMI_line <- as.factor(meta_test$PPMI_line)
  meta_test$treatment_dose <- factor(meta_test$treatment_dose, levels = c("CTL_0", treatment))
  
  dge <- DGEList(cts_test)
  dge <- calcNormFactors(dge, method = "TMM")
  
  form <- ~ treatment_dose + (1 | PPMI_line)
  
  vobjDream <- voomWithDreamWeights(dge, form, meta_test)
  fitmm <- dream(vobjDream, form, meta_test)
  fitmm <- eBayes(fitmm)
  
  res <- topTable(fitmm, coef = paste0("treatment_dose", treatment), number = Inf)
  
  return(list(res_df = res, res_obj = fitmm))
}

# --- NEW: Run the analyses in parallel using bplapply ---
# MODIFIED LINE 2: Pass the data frames to bplapply
parallel_results <- bplapply(treatments, run_dream_for_treatment, 
                             local_meta = meta_filtered, 
                             local_cts = cts_filtered)

# --- Reformat the parallel results back into your original list structures ---
res_list <- lapply(parallel_results, `[[`, "res_df")
res_obj_list <- lapply(parallel_results, `[[`, "res_obj")

names(res_list) <- treatments
names(res_obj_list) <- treatments

saveRDS(res_list, paste0("./analysis/cellculture/iMG_bulk_RNA/tmp_results/", batch, "_dream_res_df_list.rds"))
saveRDS(res_obj_list, paste0("./analysis/cellculture/iMG_bulk_RNA/tmp_results/", batch, "_dream_res_obj_list.rds"))

print(paste0("Done with batch: ", batch))
library(tidyverse)
library(RUVSeq)
library(DESeq2)
library(MAST)
library(edgeR)
library(limma)
library(nebula)
library(variancePartition)
library(ggplot2)
library(inflection)

setwd("/data/ADRD/glia_across_NDDs")

set.seed(12345)

########################################

# set necessary variables

args <- commandArgs(trailingOnly = T)
test <- args[1]

# test <- "Gerrits_TC/MTG_FTDGRN"  # dataset_region_disease
components <- unlist(strsplit(test, "_"))
dataset <- components[1]
brain_region <- components[2]
disease <- components[3]
rm(components)
rm(test)

########################################

# load in donor & cell-level metadata

donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")
donor_meta$donor <- gsub("-", ".", donor_meta$donor)
donor_meta$donor <- gsub("_", ".", donor_meta$donor)
donor_meta$group <- gsub("-", "", donor_meta$group)


celllevel_meta <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata.rds")

########################################

# load raw counts matrix (cell-level)

if (dataset == "Pineda") {
  cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Pineda_microglia_sc_cts.rds")
} else if (dataset == "AMP-PD") {
  cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/AMP-PD_microglia_sc_cts.rds")
} else if (dataset == "Mathys") {
  cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Mathys_microglia_sc_cts.rds")
} else {
  cts <- readRDS("./analysis/microglia/differential_expression/single_cell_cts_matrices/Gerrits_microglia_sc_cts.rds")
}

########################################

# filter metadata and counts to only have cells in the dataset, region, disease comparison

donor_meta <- donor_meta[donor_meta$study == dataset & donor_meta$group %in% c("HC", disease), ]

celllevel_meta$donor <- gsub("-", ".", celllevel_meta$donor)
celllevel_meta$donor <- gsub("_", ".", celllevel_meta$donor)
celllevel_meta <- celllevel_meta[rownames(celllevel_meta) %in% colnames(cts) & 
                                   celllevel_meta$region == brain_region & 
                                   celllevel_meta$donor %in% donor_meta$donor, ]


donor_meta <- donor_meta[donor_meta$donor %in% celllevel_meta$donor, ]


cts <- cts[, colnames(cts) %in% rownames(celllevel_meta)]

########################################

# pseudobulk by cluster x donor x region and filter for expressed genes

pb_cts <- aggregate(as.matrix(t(cts)), by = list(celllevel_meta$donor), FUN = sum)
colnames(pb_cts)[1] <- "donor"
pb_cts <- pb_cts %>%
  column_to_rownames(var = "donor") %>%
  t(.) %>%
  as.data.frame(.)

pb_dge <- DGEList(pb_cts)
pb_dge <- calcNormFactors(pb_dge, method = "TMM")
pb_cpm <- as.data.frame(edgeR::cpm(pb_dge))

# filter genes the same way Nebula does -- expressed in 0.5% of cells (proportion = 0.005)
# filter both single-cell and pseudobulk counts matrices

genes_to_keep <- rowSums(cts > 0) >= (ncol(cts) * 0.005)
cts_filtered <- cts[genes_to_keep, ]

pb_cts_filtered <- pb_cts[rownames(pb_cts) %in% rownames(cts_filtered), ]
pb_cpm_filtered <- pb_cpm[rownames(pb_cpm) %in% rownames(cts_filtered), ]

########################################

# run RUVseq and extract top n components of unwanted variance

rownames(donor_meta) <- NULL

RUV_meta <- donor_meta %>% 
  column_to_rownames(var = "donor")

RUV_formula <- as.formula("~ group") 

RUV_design <- model.matrix(RUV_formula, data = RUV_meta)

d_e = DGEList(pb_cts_filtered, genes = rownames(pb_cts_filtered))
d_e = calcNormFactors(d_e, method = "TMM")
d_e = estimateGLMCommonDisp(d_e, RUV_design)
d_e = estimateGLMTagwiseDisp(d_e, RUV_design)
fit1 = glmFit(d_e, RUV_design)
res1 = residuals(fit1, type = "deviance")

k_test <- ceiling(nrow(RUV_meta) / 2)
RUV_cov = RUVr(round(d_e$counts), as.character(rownames(d_e$counts)), k = k_test, res1)

RUV_cov_w <- as.data.frame(RUV_cov$W)
rownames(RUV_cov_w) <- rownames(RUV_meta)
RUV_cov_w <- RUV_cov_w %>%
  rownames_to_column(var = "donor")

########################################

# merge RUV covars w/ donor and cell-level metadata

donor_meta_merged <- donor_meta %>%
  left_join(RUV_cov_w, by = "donor")
donor_meta_merged$group <- factor(donor_meta_merged$group, levels = c("HC", disease))
donor_meta_merged$sex <- factor(donor_meta_merged$sex, levels = c("M", "F"))


celllevel_meta_merged <- celllevel_meta %>% 
  left_join(donor_meta_merged, by = "donor")
rownames(celllevel_meta_merged) <- rownames(celllevel_meta)

# make sure that order of donors in metadata matches order in counts matrix
donor_meta_merged <- donor_meta_merged %>%
  column_to_rownames(var = "donor")
donor_meta_merged <- donor_meta_merged[colnames(pb_cpm_filtered), ]
donor_meta_merged <- donor_meta_merged %>%
  rownames_to_column(var = "donor")

########################################

# run variancePartition to find inflection point for % variance explained by residuals, i.e. how many PCs of unwanted variance to regress out

varPart_meta <- donor_meta_merged %>%
  column_to_rownames(var = "donor")


# iteratively loop and run variancePartition, adding an additional term of unwanted variance every time

base_formula <- "~ group"
additional_terms <- ""
residual_variance_df <- data.frame(RUV_W = integer(), 
                                   mean_residual = numeric(), 
                                   mean_group = numeric())

for (i in 1:k_test) {
  additional_terms <- paste(additional_terms, paste0("+ W_", i))
  current_formula <- as.formula(paste(base_formula, additional_terms))
  print(current_formula)
  varPart <- fitExtractVarPartModel(pb_cpm_filtered, current_formula, varPart_meta)
  mean_residual <- mean(varPart$Residuals)
  mean_group <- mean(varPart$group)
  residual_variance_df <- rbind(residual_variance_df, data.frame(RUV_W = i, 
                                                                 mean_residual = mean_residual, 
                                                                 mean_group = mean_group))
}

# use point of maximal inflection to determine optimal number of components for regression

elbow_result <- ede(residual_variance_df$RUV_W, residual_variance_df$mean_residual, 0)
optimal_RUV_W <- elbow_result[[1]]

########################################

# create the final design formulas for DESeq and Nebula using the correct number of terms

# for testing: optimal_RUV_W = 10

terms <- paste0("W_", 1:optimal_RUV_W)

base_formula_DESeq <- "~ group"
final_formula_DEseq <- formula(paste(base_formula_DESeq, paste(terms, collapse = " + "), sep = " + "))


base_formula_nebula <- "~ group + pct.mito"
final_formula_nebula <- formula(paste(base_formula_nebula, paste(terms, collapse = " + "), sep = " + "))

########################################

# save key pieces of info from covariate selection

covar_selection_list <- list(RUV_df = RUV_cov_w, 
                             varPart_df = residual_variance_df,
                             optimal_RUV_W = optimal_RUV_W,
                             final_formula_DEseq = final_formula_DEseq,
                             final_formula_nebula = final_formula_nebula)

brain_region <- gsub("/", "_", brain_region)
saveRDS(covar_selection_list, file = paste0("./analysis/microglia/differential_expression/covar_selection/", dataset, "_", brain_region, "_", disease, "_covar_selection.rds"))

########################################
########################################
########################################

# run DESeq2

dds <- DESeqDataSetFromMatrix(countData = pb_cts_filtered,
                              colData = donor_meta_merged, 
                              design = final_formula_DEseq)

dds$group <- relevel(dds$group, ref = "HC")

dds <- DESeq(dds)
DESeq_results <- as.data.frame(results(dds, name = paste0("group_", disease, "_vs_HC")))
DESeq_results <- na.omit(DESeq_results)
DESeq_results <- DESeq_results %>%
  rownames_to_column(var = "gene")

saveRDS(DESeq_results, file = paste0("./analysis/microglia/differential_expression/DESeq_results/", dataset, "_", brain_region, "_", disease, "_DESeq_results.rds"))

########################################
########################################
########################################

# run nebula

nebula_donor_meta <- donor_meta_merged %>%
  column_to_rownames(var = "donor")

nebula_donor_idents <- celllevel_meta_merged$donor


nebula <- list(count = cts_filtered,
               pred = celllevel_meta_merged, 
               sid = nebula_donor_idents,
               offset = celllevel_meta_merged$nCount_RNA)

nebula_design <- model.matrix(final_formula_nebula, data = nebula$pred)

# if the cells are already grouped, group_cell() won't run, and need to use nebula object above
# this takes the text that is printed ("the cells are already grouped") and writes an if/else statement with it
# if text is printed when running group_cell(), runs nebula() with nebula object above
# if text is not printed, runs group_cell() and runs nebula() with the output
nebula_group_output <- capture.output(result <- group_cell(count = nebula$count, id = nebula$sid, pred = nebula_design, offset = nebula$pred$nCount_RNA))

if (any(grepl("The cells are already grouped", nebula_group_output))) {
  nebula_results <- nebula(nebula$count, nebula$sid, pred = nebula_design, 
                           offset = nebula$offset, ncore = 1)
} else {
  nebula_grouped = group_cell(count = nebula$count, id = nebula$sid, pred = nebula_design, offset = nebula$offset)
  nebula_results <- nebula(nebula_grouped$count, nebula_grouped$id, pred = nebula_design, 
                           offset = nebula_grouped$offset, ncore = 1)
}

saveRDS(nebula_results, file = paste0("./analysis/microglia/differential_expression/nebula_results/", dataset, "_", brain_region, "_", disease, "_nebula_results.rds"))

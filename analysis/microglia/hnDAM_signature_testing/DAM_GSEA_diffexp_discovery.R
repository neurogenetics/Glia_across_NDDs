library(tidyverse)
library(fgsea)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
path <- list(DAM_genes)
names(path) <- "DAM"

diffexp_tests <- read.delim("./analysis/diffexp_tests.txt", header = F)
diffexp_tests <- gsub("/", "_", diffexp_tests$V1)
#test = diffexp_tests[2]

##################################################
##################################################
##################################################

gsea_summary <- data.frame()
rankings <- list()

for (test in diffexp_tests){
  nebula <- readRDS(paste0("./analysis/microglia/differential_expression/nebula_results/", test, "_nebula_results.rds"))
  disease <- sub(".*_(.*)$", "\\1", test)
  
  summary <- nebula$summary
  summary <- summary[c(paste0("logFC_group", disease), paste0("p_group", disease), paste0("se_group", disease), "gene")]
  summary$convergence <- nebula$convergence
  summary <- summary[summary$convergence == 1, ]
  summary <- summary[summary[[3]] < 10, ]
  summary$rank <- sign(summary[[1]]) * -log10(summary[[2]])
  summary <- summary %>%
    arrange(desc(rank))
  
  ranking <- summary$rank
  names(ranking) <- summary$gene
  
  gseares <- fgsea(pathways = path, stats = ranking)
  gseares$test <- test
  gseares <- as.data.frame(gseares)
  
  rankings[[test]] <- ranking

  gsea_summary <- rbind(gseares, gsea_summary)
}

gsea_summary$padj <- p.adjust(gsea_summary$pval, method = "BH")

saveRDS(gsea_summary, file = "./analysis/microglia/DAM_signature_discovery/DAM_diffexp_GSEA_summary.rds")
saveRDS(rankings, file = "./analysis/microglia/DAM_signature_discovery/DAM_diffexp_GSEA_gene_rankings.rds")

library(tidyverse)
library(ComplexUpset)
library(eulerr)
library(svglite)
library(fgsea)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

hnDAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")


micro_signatures <- read.csv("./analysis/microglia/misc_files/micro_DAM_signatures_to_compare.csv", header = T)
micro_signatures_list <- lapply(micro_signatures, function(x) x[x != ""])
names(micro_signatures_list) <- names(micro_signatures)

##################################################

micro_signatures_list[["hnDAM"]] <- hnDAM_genes

micro_signatures_list <- micro_signatures_list[c("hnDAM", "Sun_2023_MG4_human", "Marshe_2025_GPNMB.hi_human", 
                                                 "Silvin_2022_human_DAMs_human" ,"Keren.Shaul_2017_DAM_mouse", 
                                                 "Krasemann_2017_MGnD_mouse")]

names(micro_signatures_list) <- c("hnDAM", "Sun_2023_MG4", "Marshe_2025_GPNMB-hi", 
                                  "Silvin_2022_human_DAMs" ,"Keren-Shaul_2017_DAM_5xFAD", 
                                  "Krasemann_2017_MGnD_APP-PS1")

all_items <- unique(unlist(micro_signatures_list))
binary_matrix <- sapply(micro_signatures_list, function(x) all_items %in% x)
binary_df <- as.data.frame(binary_matrix)
binary_df <- cbind(item = all_items, binary_df)
binary_df <- binary_df %>%
  column_to_rownames(var = "item")

binary_df[] <- lapply(binary_df, as.integer)



ComplexUpset::upset(binary_df, 
                    rev(colnames(binary_df)[!grepl("DAM2", colnames(binary_df))]), 
                    mode = "inclusive_intersection", 
                    # intersections = list(c("hnDAM", "Sun_2023_MG4_human"),
                    #                      c("hnDAM", "Marshe_2025_GPNMB-hi_human"),
                    #                      c("hnDAM", "Silvin_2022_human_DAMs_human"),
                    #                      c("hnDAM", "Keren-Shaul_2017_DAM_mouse"),
                    #                      c("hnDAM", "Krasemann_2017_MGnD_mouse")), 
                    set_sizes = FALSE,
                    sort_sets = FALSE)

##################################################
##################################################
##################################################

fit <- venn(micro_signatures_list[c(4, 2, 3, 1)])

png("./analysis/microglia/plots_final/DAM_overlap_human_venn.png", width = 3, height = 3, units = "in", res = 600)
plot(fit, fills = list(fill = c("#de76d3", "#910c19", "#7a449c", "yellow"), alpha = 0.6),
     labels = F, edges = T, quantities = T)
dev.off()



fit2 <- venn(micro_signatures_list[c(6, 5, 1)])

png("./analysis/microglia/plots_final/DAM_overlap_mouse_venn.png", width = 3, height = 3, units = "in", res = 600)
plot(fit2, fills = list(fill = c("#910c19", "#7a449c", "yellow"), alpha = 0.6),
     labels = F, edges = T, quantities = T)
dev.off()

##################################################
##################################################
##################################################

# check intersecting genes

y <- sort(intersect(micro_signatures_list$hnDAM, 
          intersect(micro_signatures_list$Sun_2023_MG4, micro_signatures_list$`Marshe_2025_GPNMB-hi`)))

x <- sort(intersect(micro_signatures_list$hnDAM, micro_signatures_list$Sun_2023_MG4_human))

z <- x[!x %in% y]
z


sort(intersect(micro_signatures_list$hnDAM, 
               intersect(micro_signatures_list$Sun_2023_MG4, 
                         intersect(micro_signatures_list$`Marshe_2025_GPNMB-hi`, micro_signatures_list$Silvin_2022_human_DAMs))))


sort(intersect(micro_signatures_list$hnDAM, 
                    union(micro_signatures_list$`Keren-Shaul_2017_DAM_5xFAD`, micro_signatures_list$`Krasemann_2017_MGnD_APP-PS1`)))

##################################################
##################################################
##################################################

GEPs <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_GEPs_zscore_df.rds")

ranking <- GEPs$avg_z
names(ranking) <- GEPs$gene


gseares <- fgsea(pathways = micro_signatures_list, stats = ranking, eps = 0)

plotEnrichment(pathway = micro_signatures_list$Sun_2023_MG4_human, stats = ranking)

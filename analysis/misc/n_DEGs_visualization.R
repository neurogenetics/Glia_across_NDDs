library(tidyverse)
library(ggplot2)
library(patchwork)
library(svglite)

setwd("/data/ADRD/glia_across_NDDs")

########################################
########################################
########################################

astro_res_list <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_diffexp_results_for_mashr.rds")
micro_res_list <- readRDS("./analysis/microglia/differential_expression/microglia_diffexp_results_for_mashr.rds")
oligo_res_list <- readRDS("./analysis/oligodendrocytes/differential_expression/oligodendrocytes_diffexp_results_for_mashr.rds")

source("./code_organized/functions/rename_datasets.R")

########################################
########################################
########################################

# padjust nebula results

n_DEGs_df <- data.frame()

for (test in names(astro_res_list)){
  res = astro_res_list[[test]]
  res$padj_nebula <- p.adjust(res[[4]])
  astro_res_list[[test]] <- res
  
  df <- data.frame(celltype = "Astrocytes",
                   test = test,
                   up_nebula = sum(res$padj_nebula < 0.05 & res[[2]] > 0),
                   down_nebula = sum(res$padj_nebula < 0.05 & res[[2]] < 0),
                   up_DESeq = sum(res$padj < 0.05 & res[[6]] > 0),
                   down_DESeq = sum(res$padj < 0.05 & res[[6]] < 0),
                   
                   discordant = sum(res$padj_nebula < 0.05 & res$padj < 0.05 &
                                      !sign(res[[2]]) == sign(res[[6]])))
  n_DEGs_df <- rbind(n_DEGs_df, df)
}


for (test in names(micro_res_list)){
  res = micro_res_list[[test]]
  res$padj_nebula <- p.adjust(res[[4]])
  micro_res_list[[test]] <- res
  
  df <- data.frame(celltype = "Microglia",
                   test = test,
                   up_nebula = sum(res$padj_nebula < 0.05 & res[[2]] > 0),
                   down_nebula = sum(res$padj_nebula < 0.05 & res[[2]] < 0),
                   up_DESeq = sum(res$padj < 0.05 & res[[6]] > 0),
                   down_DESeq = sum(res$padj < 0.05 & res[[6]] < 0),
                   
                   discordant = sum(res$padj_nebula < 0.05 & res$padj < 0.05 &
                                      !sign(res[[2]]) == sign(res[[6]])))
  n_DEGs_df <- rbind(n_DEGs_df, df)
}


for (test in names(oligo_res_list)){
  res = oligo_res_list[[test]]
  res$padj_nebula <- p.adjust(res[[4]])
  oligo_res_list[[test]] <- res
  
  df <- data.frame(celltype = "Oligodendrocytes",
                   test = test,
                   up_nebula = sum(res$padj_nebula < 0.05 & res[[2]] > 0),
                   down_nebula = sum(res$padj_nebula < 0.05 & res[[2]] < 0),
                   up_DESeq = sum(res$padj < 0.05 & res[[6]] > 0),
                   down_DESeq = sum(res$padj < 0.05 & res[[6]] < 0),
                   
                   discordant = sum(res$padj_nebula < 0.05 & res$padj < 0.05 &
                                      !sign(res[[2]]) == sign(res[[6]])))
  n_DEGs_df <- rbind(n_DEGs_df, df)
}


n_DEGs_df$test <- rename_datasets(n_DEGs_df$test)

########################################
########################################
########################################

# get # cells/donors for each comparison

micro_celllevel_meta <- readRDS("./analysis/microglia/all_microglia_celllevel_metadata_ANNOTATED.rds")
astro_celllevel_meta <- readRDS("./analysis/astrocytes/all_astros_metadata_ANNOTATED.rds")
oligos_celllevel_meta <- readRDS("./analysis/oligodendrocytes/oligos_celllevel_metadata_ANNOTATED.rds")
oligos_celllevel_meta$region <- ifelse(oligos_celllevel_meta$region == "HC", "HIP", oligos_celllevel_meta$region)


donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")
donor_meta <- readRDS("./metadata/cleaned/meta_merged.rds")
donor_meta$donor <- gsub("-", ".", donor_meta$donor)
donor_meta$donor <- gsub("_", ".", donor_meta$donor)
donor_meta$group <- gsub("-", "", donor_meta$group)


tests <- read.delim("./analysis/diffexp_tests.txt", header = F)
tests <- tests$V1
tests_oligos <- read.delim("./analysis/diffexp_tests_oligos.txt", header = F)
tests_oligos <- tests_oligos$V1


summary <- data.frame()

for (test in tests){
  components <- unlist(strsplit(test, "_"))
  dataset <- components[1]
  brain_region <- components[2]
  disease <- components[3]
  
  donor_meta_filtered <- donor_meta[donor_meta$study == dataset & donor_meta$group %in% c("HC", disease), ]
  
  celllevel_meta <- micro_celllevel_meta
  celllevel_meta$donor <- gsub("-", ".", celllevel_meta$donor)
  celllevel_meta$donor <- gsub("_", ".", celllevel_meta$donor)
  celllevel_meta <- celllevel_meta[celllevel_meta$dataset == dataset & 
                                     celllevel_meta$region == brain_region & 
                                     celllevel_meta$donor %in% donor_meta_filtered$donor, ]
  celllevel_meta <- celllevel_meta %>%
    left_join(donor_meta_filtered, by = "donor")
  
  df <- data.frame(celltype = "Microglia",
                   test = test,
                   n_cells_HC = sum(celllevel_meta$group == "HC"),
                   n_cells_disease = sum(celllevel_meta$group == disease),
                   n_donors_HC = length(unique(celllevel_meta$donor[celllevel_meta$group == "HC"])),
                   n_donors_disease = length(unique(celllevel_meta$donor[celllevel_meta$group == disease])),
                   med_UMIs = median(celllevel_meta$nCount_RNA))
  
  summary <- rbind(summary, df)
}

for (test in tests){
  components <- unlist(strsplit(test, "_"))
  dataset <- components[1]
  brain_region <- components[2]
  disease <- components[3]
  
  donor_meta_filtered <- donor_meta[donor_meta$study == dataset & donor_meta$group %in% c("HC", disease), ]
  
  celllevel_meta <- astro_celllevel_meta
  celllevel_meta$donor <- gsub("-", ".", celllevel_meta$donor)
  celllevel_meta$donor <- gsub("_", ".", celllevel_meta$donor)
  celllevel_meta <- celllevel_meta[celllevel_meta$dataset == dataset & 
                                     celllevel_meta$region == brain_region & 
                                     celllevel_meta$donor %in% donor_meta_filtered$donor, ]
  celllevel_meta <- celllevel_meta %>%
    left_join(donor_meta_filtered, by = "donor")
  
  df <- data.frame(celltype = "Astrocytes",
                   test = test,
                   n_cells_HC = sum(celllevel_meta$group == "HC"),
                   n_cells_disease = sum(celllevel_meta$group == disease),
                   n_donors_HC = length(unique(celllevel_meta$donor[celllevel_meta$group == "HC"])),
                   n_donors_disease = length(unique(celllevel_meta$donor[celllevel_meta$group == disease])),
                   med_UMIs = median(celllevel_meta$nCount_RNA))
  
  summary <- rbind(summary, df)
}

for (test in tests_oligos){
  components <- unlist(strsplit(test, "_"))
  dataset <- components[1]
  brain_region <- components[2]
  disease <- components[3]
  
  donor_meta_filtered <- donor_meta[donor_meta$study == dataset & donor_meta$group %in% c("HC", disease), ]
  
  celllevel_meta <- oligos_celllevel_meta
  celllevel_meta$donor <- gsub("-", ".", celllevel_meta$donor)
  celllevel_meta$donor <- gsub("_", ".", celllevel_meta$donor)
  celllevel_meta <- celllevel_meta[celllevel_meta$dataset == dataset & 
                                     celllevel_meta$region == brain_region & 
                                     celllevel_meta$donor %in% donor_meta_filtered$donor, ]
  celllevel_meta <- celllevel_meta %>%
    left_join(donor_meta_filtered, by = "donor")
  
  df <- data.frame(celltype = "Oligodendrocytes",
                   test = test,
                   n_cells_HC = sum(celllevel_meta$group == "HC"),
                   n_cells_disease = sum(celllevel_meta$group == disease),
                   n_donors_HC = length(unique(celllevel_meta$donor[celllevel_meta$group == "HC"])),
                   n_donors_disease = length(unique(celllevel_meta$donor[celllevel_meta$group == disease])),
                   med_UMIs = median(celllevel_meta$nCount_RNA))
  
  summary <- rbind(summary, df)
}


summary$test <- rename_datasets(summary$test)

saveRDS(summary, file = "./analysis/misc/cell_donor_counts_for_DE_summary_fig.rds")

########################################
########################################
########################################

# make plots

df_nebula <- n_DEGs_df %>%
  dplyr::select(celltype, test, up_nebula, down_nebula) %>%
  pivot_longer(cols = c(up_nebula, down_nebula), names_to = "direction", values_to = "count") %>%
  mutate(count = ifelse(direction == "down_nebula", -count, count))


nebula_astros <- df_nebula[df_nebula$celltype == "Astrocytes", ]
nebula_micro <- df_nebula[df_nebula$celltype == "Microglia", ]
nebula_oligos <- df_nebula[df_nebula$celltype == "Oligodendrocytes", ]


p1 <- ggplot(nebula_astros, aes(x = test, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +
  labs(title = "nebula") +
  scale_fill_manual(values = c("up_nebula" = "firebrick", "down_nebula" = "steelblue")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

p2 <- ggplot(nebula_astros, aes(x = test, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +
  labs(title = "nebula") +
  scale_fill_manual(values = c("up_nebula" = "firebrick", "down_nebula" = "steelblue")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

p3 <- ggplot(nebula_oligos, aes(x = test, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +
  labs(title = "nebula") +
  scale_fill_manual(values = c("up_nebula" = "firebrick", "down_nebula" = "steelblue")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

########################################

df_DESeq <- n_DEGs_df %>%
  dplyr::select(celltype, test, up_DESeq, down_DESeq) %>%
  pivot_longer(cols = c(up_DESeq, down_DESeq), names_to = "direction", values_to = "count") %>%
  mutate(count = ifelse(direction == "down_DESeq", -count, count))

DESeq_astros <- df_DESeq[df_DESeq$celltype == "Astrocytes", ]
DESeq_micro <- df_DESeq[df_DESeq$celltype == "Microglia", ]
DESeq_oligos <- df_DESeq[df_DESeq$celltype == "Oligodendrocytes", ]


p4 <- ggplot(DESeq_astros, aes(x = test, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +
  labs(title = "DESeq") +
  scale_fill_manual(values = c("up_DESeq" = "firebrick", "down_DESeq" = "steelblue")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

p5 <- ggplot(DESeq_astros, aes(x = test, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +
  labs(title = "DESeq") +
  scale_fill_manual(values = c("up_DESeq" = "firebrick", "down_DESeq" = "steelblue")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

p6 <- ggplot(DESeq_oligos, aes(x = test, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +
  labs(title = "DESeq") +
  scale_fill_manual(values = c("up_DESeq" = "firebrick", "down_DESeq" = "steelblue")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

########################################

# number of cells / donors
summary <- readRDS("./analysis/misc/cell_donor_counts_for_DE_summary_fig.rds")

df_ncells <- summary %>%
  dplyr::select(celltype, test, n_cells_HC, n_cells_disease) %>%
  pivot_longer(cols = c(n_cells_HC, n_cells_disease), names_to = "group", values_to = "count")

ncells_astros <- df_ncells[df_ncells$celltype == "Astrocytes", ]
ncells_micro <- df_ncells[df_ncells$celltype == "Microglia", ]
ncells_oligos <- df_ncells[df_ncells$celltype == "Oligodendrocytes", ]


p7 <- ggplot(ncells_astros, aes(x = test, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "# cells") +
  scale_fill_manual(values = c("n_cells_disease" = "#ac54bf", "n_cells_HC" = "grey75")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))


p8 <- ggplot(ncells_micro, aes(x = test, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "# cells") +
  scale_fill_manual(values = c("n_cells_disease" = "#ac54bf", "n_cells_HC" = "grey75")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))


p9 <- ggplot(ncells_oligos, aes(x = test, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "# cells") +
  scale_fill_manual(values = c("n_cells_disease" = "#ac54bf", "n_cells_HC" = "grey75")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

########################################

# number of donors

df_ndonors <- summary %>%
  dplyr::select(celltype, test, n_donors_HC, n_donors_disease) %>%
  pivot_longer(cols = c(n_donors_HC, n_donors_disease), names_to = "group", values_to = "count")

ndonors_astros <- df_ndonors[df_ndonors$celltype == "Astrocytes", ]
ndonors_micro <- df_ndonors[df_ndonors$celltype == "Microglia", ]
ndonors_oligos <- df_ndonors[df_ndonors$celltype == "Oligodendrocytes", ]


p10 <- ggplot(ndonors_astros, aes(x = test, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "# donors") +
  scale_fill_manual(values = c("n_donors_disease" = "#ac54bf", "n_donors_HC" = "grey75")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))


p11 <- ggplot(ndonors_micro, aes(x = test, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "# donors") +
  scale_fill_manual(values = c("n_donors_disease" = "#ac54bf", "n_donors_HC" = "grey75")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))


p12 <- ggplot(ndonors_oligos, aes(x = test, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "# donors") +
  scale_fill_manual(values = c("n_donors_disease" = "#ac54bf", "n_donors_HC" = "grey75")) +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

########################################

# median UMIs

df_UMIs <- summary %>%
  dplyr::select(celltype, test, med_UMIs) 


UMIs_astros <- df_UMIs[df_UMIs$celltype == "Astrocytes", ]
UMIs_micro <- df_UMIs[df_UMIs$celltype == "Microglia", ]
UMIs_oligos <- df_UMIs[df_UMIs$celltype == "Oligodendrocytes", ]


p13 <- ggplot(UMIs_astros, aes(x = test, y = med_UMIs, fill = med_UMIs)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "median UMIs") +
  scale_fill_gradient2() +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))


p14 <- ggplot(UMIs_micro, aes(x = test, y = med_UMIs, fill = med_UMIs)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "median UMIs") +
  scale_fill_gradient2() +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))



p15 <- ggplot(UMIs_oligos, aes(x = test, y = med_UMIs, fill = med_UMIs)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "median UMIs") +
  scale_fill_gradient2() +
  theme_bw() + 
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "None",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))


########################################

p16 <- p1 + p4 + p7 + p10 + p13 + plot_layout(nrow = 1)
p17 <- p2 + p5 + p8 + p11 + p14 + plot_layout(nrow = 1)
p18 <- p3 + p6 + p9 + p12 + p15 + plot_layout(nrow = 1)

p19 <- p16 / p17 / p18

svglite("./analysis/misc/DE_cells_summary_fig.svg", width = 14.5, height = 14.5)
plot(p19)
dev.off()

ggsave(plot = p19, filename = "./analysis/misc/DE_cells_summary_fig.png", width = 14.5, height = 14.5, units = "in", dpi = 600)

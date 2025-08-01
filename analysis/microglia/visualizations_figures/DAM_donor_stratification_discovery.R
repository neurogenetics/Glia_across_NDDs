library(tidyverse)
library(ggplot2)
library(svglite)
library(patchwork)
library(fgsea)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

source("./code_organized/functions/plotting_colors.R")

##################################################
##################################################
##################################################

# summary fig for DAM replication -- GSEA + PCA

PCA_summary <- readRDS("./analysis/microglia/DAM_signature_discovery/PCA_GLM_summary.rds")
gsea_summary <- readRDS("./analysis/microglia/DAM_signature_discovery/DAM_diffexp_GSEA_summary.rds")

##################################################

PCA_summary$test <- c("Mathys_2024_EC_AD", "Mathys_2024_EC_CR", "Mathys_2024_HIP_AD", "Mathys_2024_HIP_CR", 
                      "Mathys_2024_MTG_AD", "Mathys_2024_MTG_CR", "Mathys_2024_AnG_AD", "Mathys_2024_AnG_CR", 
                      "Mathys_2024_TH_AD", "Mathys_2024_TH_CR", "Mathys_2024_PFC_AD", "Mathys_2024_PFC_CR", 
                      "NM_2024_DMV_PD", "NM_2024_GPi_PD", "NM_2024_M1_PD", "NM_2024_PFC_PD", "NM_2024_V1_PD",
                      "Pineda_2024_M1_ALS", "Pineda_2024_M1_C9ALS", "Pineda_2024_M1_FTLD", "Pineda_2024_M1_C9FTLD", 
                      "Pineda_2024_PFC_ALS", "Pineda_2024_PFC_C9ALS", "Pineda_2024_PFC_FTLD", "Pineda_2024_PFC_C9FTLD", 
                      "Gerrits_2022_TC_FTD-GRN", "Gerrits_2022_OC_FTD-GRN", "Gerrits_2022_FC_FTD-GRN")

gsea_summary$test <- rev(PCA_summary$test)

##################################################

gsea_summary <- gsea_summary %>%
  arrange(desc(NES))
gsea_summary$test <- factor(gsea_summary$test, levels = rev(gsea_summary$test))

gsea_summary$`FDR < 0.05` <- ifelse(gsea_summary$padj < 0.05, T, F)

gsea_summary$`Disease context` <- ifelse(grepl("Mathys", gsea_summary$test), "AD spectrum",
                                         ifelse(grepl("Pineda|Gerrits", gsea_summary$test), "ALS-FTD spectrum",
                                                ifelse(grepl("NM", gsea_summary$test), "PD-DLB spectrum", "")))


p1 <- ggplot(gsea_summary, aes(x = NES, y = test)) +
  geom_segment(aes(x = 0, xend = NES, y = test, yend = test, color = `FDR < 0.05`),
               size = 1, show.legend = TRUE) +
  geom_point(aes(fill = `Disease context`, shape = `Measure type`, color = `FDR < 0.05`),
             size = 5, stroke = 1.2, show.legend = TRUE, shape = 22) +
  scale_fill_manual(values = disease_colors) +
  scale_color_manual(name = "FDR < 0.05", values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  guides(color = guide_legend(override.aes = list(shape = 21, fill = "white", size = 5, stroke = 1.2)),
         fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 5)),
         shape = guide_legend()) +
  theme_bw() +
  labs(x = "NES") +
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

#########################

PCA_summary$Estimate <- abs(PCA_summary$Estimate)

PCA_summary$test <- factor(PCA_summary$test, levels = rev(gsea_summary$test))

PCA_summary$test <- factor(PCA_summary$test, levels = rev(gsea_summary$test))

PCA_summary$`FDR < 0.05` <- ifelse(PCA_summary$padj < 0.05, T, F)

PCA_summary$`Disease context` <- ifelse(grepl("Mathys", PCA_summary$test), "AD spectrum", 
                                        ifelse(grepl("Pineda|Gerrits", PCA_summary$test), "ALS-FTD spectrum", 
                                               ifelse(grepl("NM", PCA_summary$test), "PD-DLB spectrum", "")))

PCA_summary$`Measure type` <- "Clinical"


p2 <- ggplot(PCA_summary, aes(x = Estimate, y = test)) +
  geom_segment(aes(x = 0, xend = Estimate, y = test, yend = test, color = `FDR < 0.05`),
               size = 1, show.legend = TRUE) +
  geom_point(aes(fill = `Disease context`, shape = `Measure type`, color = `FDR < 0.05`),
             size = 5, stroke = 1.2, show.legend = TRUE) +
  scale_shape_manual(values = c("Pathological" = 21, "Clinical" = 22)) +
  scale_fill_manual(values = disease_colors) +
  scale_color_manual(name = "FDR < 0.05", values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
  guides(color = guide_legend(override.aes = list(shape = 21, fill = "white", size = 5, stroke = 1.2)),
         fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 5)),
         shape = guide_legend()) +
  theme_bw() +
  labs(x = "| GLM estimate |") +
  theme(legend.position = "right",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

#########################

svglite(filename = "./analysis/microglia/plots_final/DAM_discovery_summary_plot.svg", width = 9.5, height = 8.5)
plot(p1 + p2)
dev.off()

##################################################
##################################################
##################################################

# plots_final of individual GSEA & PCA

# Mathys 2024: HIP / AD
# Pineda 2024: M1 / FTLD
# NM 2024: PFC
# Gerrits 2022: TC

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
pathway <- list(DAM_genes)
names(pathway) <- "DAM"

PCA_df_list <- readRDS("./analysis/microglia/DAM_signature_discovery/PCA_df_list.rds")
GSEA_ranking_list <- readRDS("./analysis/microglia/DAM_signature_discovery/DAM_diffexp_GSEA_gene_rankings.rds")
variance_explained_list <- readRDS("./analysis/microglia/DAM_signature_discovery/variance_explained_list.rds")

##################################################
##################################################
##################################################

# Mathys 2024 - HIP

##################################################

Mathys_PCA <- PCA_df_list[["Mathys_HIP_AD"]]

p3 <- ggplot(Mathys_PCA, aes(x = group, y = PC1, fill = group)) +
  geom_boxplot() + 
  geom_point() + 
  scale_fill_manual(values = disease_colors) +
  theme_bw() +
  labs(x = "Diagnosis",
       y = "PC1 (30.80%)") + 
  scale_y_continuous(limits = c(-9, 9)) +
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

svglite(filename = "./analysis/microglia/plots_final/Mathys_2024_HIP_AD_PCA.svg", width = 6, height = 3.2)
plot(p3)
dev.off()

##################################################

ranking_Mathys_HIP <- GSEA_ranking_list[["Mathys_HIP_AD"]]

gsea_Mathys_HIP <- fgsea(pathways = pathway, stats = ranking_Mathys_HIP, scoreType = "std")

pd <- plotEnrichmentData(
  pathway = pathway$DAM,,
  stats = ranking_Mathys_HIP)
p4 <- with(pd,
           ggplot(data=curve) +
             geom_line(aes(x=rank, y=ES), color="green") +
             geom_segment(data=ticks, mapping=aes(x=rank, y=-spreadES/16, xend=rank, yend=spreadES/16), size=0.2) +
             geom_hline(yintercept=posES, colour="red", linetype="dashed") +
             geom_hline(yintercept=negES, colour="red", linetype="dashed") +
             geom_hline(yintercept=0, colour="black") +
             theme_bw() +
             theme(legend.position = "none",
                   legend.title = element_text(size = 14),
                   legend.text = element_text(size = 12),  
                   text = element_text(family = "Arial"),
                   axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12)) + 
             labs(x="Rank - AD vs. HC", y="Enrichment score"))

svglite(filename = "./analysis/microglia/plots_final/Mathys_2024_HIP_AD_GSEA.svg", width = 5.8, height = 3.2)
plot(p4)
dev.off()

##################################################
##################################################
##################################################

# Pineda 2024 - M1

##################################################

Pineda_PCA <- PCA_df_list[["Pineda_M1_FTLD"]]

p5 <- ggplot(Pineda_PCA, aes(x = group, y = PC1, fill = group)) +
  geom_boxplot() + 
  geom_point() + 
  scale_fill_manual(values = disease_colors) +
  theme_bw() +
  labs(x = "Diagnosis",
       y = "PC1 (33.05%%)") + 
  scale_y_continuous(limits = c(NA, 9.2)) +
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

svglite(filename = "./analysis/microglia/plots_final/Pineda_2024_M1_FTLD_PCA.svg", width = 6, height = 3.2)
plot(p5)
dev.off()

##################################################

ranking_Pineda_M1 <- GSEA_ranking_list[["Pineda_M1_FTLD"]]

gsea_Pineda_M1 <- fgsea(pathways = pathway, stats = ranking_Pineda_M1, scoreType = "std")

pd <- plotEnrichmentData(
  pathway = pathway$DAM,,
  stats = ranking_Pineda_M1)
p6 <- with(pd,
           ggplot(data=curve) +
             geom_line(aes(x=rank, y=ES), color="green") +
             geom_segment(data=ticks, mapping=aes(x=rank, y=-spreadES/16, xend=rank, yend=spreadES/16), size=0.2) +
             geom_hline(yintercept=posES, colour="red", linetype="dashed") +
             geom_hline(yintercept=negES, colour="red", linetype="dashed") +
             geom_hline(yintercept=0, colour="black") +
             theme_bw() +
             theme(legend.position = "none",
                   legend.title = element_text(size = 14),
                   legend.text = element_text(size = 12),  
                   text = element_text(family = "Arial"),
                   axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12)) + 
             labs(x="Rank - FTLD vs. HC", y="Enrichment score"))

svglite(filename = "./analysis/microglia/plots_final/Pineda_2024_M1_FTLD_GSEA.svg", width = 5.8, height = 3.2)
plot(p6)
dev.off()

##################################################
##################################################
##################################################

# Gerrits 2022 - TC

##################################################

Gerrits_PCA <- PCA_df_list[["Gerrits_TC/MTG_FTDGRN"]]
Gerrits_PCA$group <- ifelse(Gerrits_PCA$group == "HC", "HC", "FTD-GRN")
Gerrits_PCA$group <- factor(Gerrits_PCA$group, levels = c("HC", "FTD-GRN"))

p7 <- ggplot(Gerrits_PCA, aes(x = group, y = PC1, fill = group)) +
  geom_boxplot() + 
  geom_point() + 
  scale_fill_manual(values = disease_colors) +
  theme_bw() +
  labs(x = "Diagnosis",
       y = "PC1 (26.15%)") + 
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

svglite(filename = "./analysis/microglia/plots_final/Gerrits_2022_TC_PCA.svg", width = 6, height = 3.2)
plot(p7)
dev.off()

##################################################

ranking_Gerrits_TC <- GSEA_ranking_list[["Gerrits_TC_MTG_FTDGRN"]]

gsea_Gerrits_TC <- fgsea(pathways = pathway, stats = ranking_Gerrits_TC, scoreType = "std")

pd <- plotEnrichmentData(
  pathway = pathway$DAM,,
  stats = ranking_Gerrits_TC)
p8 <- with(pd,
           ggplot(data=curve) +
             geom_line(aes(x=rank, y=ES), color="green") +
             geom_segment(data=ticks, mapping=aes(x=rank, y=-spreadES/16, xend=rank, yend=spreadES/16), size=0.2) +
             geom_hline(yintercept=posES, colour="red", linetype="dashed") +
             geom_hline(yintercept=negES, colour="red", linetype="dashed") +
             geom_hline(yintercept=0, colour="black") +
             theme_bw() +
             theme(legend.position = "none",
                   legend.title = element_text(size = 14),
                   legend.text = element_text(size = 12),  
                   text = element_text(family = "Arial"),
                   axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12)) + 
             labs(x="Rank - FTD-GRN vs. HC", y="Enrichment score"))

svglite(filename = "./analysis/microglia/plots_final/Gerrits_2022_TC_GSEA.svg", width = 5.8, height = 3.2)
plot(p8)
dev.off()

##################################################
##################################################
##################################################

# NM 2024 - DMV

##################################################

NM_PCA <- PCA_df_list[["AMP-PD_DMV_PD"]]

p9 <- ggplot(NM_PCA, aes(x = group, y = PC1, fill = group)) +
  geom_boxplot() + 
  geom_point() + 
  scale_fill_manual(values = disease_colors) +
  theme_bw() +
  labs(x = "Diagnosis",
       y = "PC1 (15.73%)") + 
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

svglite(filename = "./analysis/microglia/plots_final/NM_2024_DMV_PCA.svg", width = 6, height = 3.2)
plot(p9)
dev.off()

##################################################

ranking_NM_DMV <- GSEA_ranking_list[["AMP-PD_DMV_PD"]]

gsea_NM_DMV <- fgsea(pathways = pathway, stats = ranking_NM_DMV, scoreType = "std")

pd <- plotEnrichmentData(
  pathway = pathway$DAM,,
  stats = ranking_NM_DMV)
p10 <- with(pd,
           ggplot(data=curve) +
             geom_line(aes(x=rank, y=ES), color="green") +
             geom_segment(data=ticks, mapping=aes(x=rank, y=-spreadES/16, xend=rank, yend=spreadES/16), size=0.2) +
             geom_hline(yintercept=posES, colour="red", linetype="dashed") +
             geom_hline(yintercept=negES, colour="red", linetype="dashed") +
             geom_hline(yintercept=0, colour="black") +
             theme_bw() +
             theme(legend.position = "none",
                   legend.title = element_text(size = 14),
                   legend.text = element_text(size = 12),  
                   text = element_text(family = "Arial"),
                   axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12)) + 
             labs(x="Rank - PD vs. HC", y="Enrichment score"))

svglite(filename = "./analysis/microglia/plots_final/NM_2024_DMV_GSEA.svg", width = 5.8, height = 3.2)
plot(p10)
dev.off()

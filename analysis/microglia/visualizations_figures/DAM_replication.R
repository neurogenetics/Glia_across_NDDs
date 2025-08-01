library(tidyverse)
library(ggplot2)
library(svglite)
library(patchwork)
library(fgsea)
library(DESeq2)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

source("./code_organized/functions/plotting_colors.R")

##################################################
##################################################
##################################################

# summary fig for DAM replication -- GSEA + PCA

gsea_summary <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_DAM_GSEA_summary.rds")
PCA_summary <- readRDS("./analysis/microglia/DAM_signature_replication/GLM_summary_merged_FDRcorr.rds")

#########################

gsea_summary <- gsea_summary %>%
  arrange(desc(NES))
gsea_summary$test <- factor(gsea_summary$test, levels = rev(gsea_summary$test))

gsea_summary$`FDR < 0.05` <- ifelse(gsea_summary$padj < 0.05, T, F)

gsea_summary$`Disease context` <- ifelse(grepl("CERAD|Braak|AD_dementia|INS_AD|MCI", gsea_summary$test), "AD spectrum", 
                                         ifelse(grepl("Macnair", gsea_summary$test), "Multiple sclerosis",
                                                ifelse(grepl("INS_PSP|INS_PiD", gsea_summary$test), "Primary tauopathy",
                                                       ifelse(grepl("ALS", gsea_summary$test), "ALS-FTD spectrum", 
                                                              ifelse(grepl("PD", gsea_summary$test), "PD-DLB spectrum", "")))))

gsea_summary$`Measure type` <- ifelse(grepl("CERAD|Braak|Rexach|Macnair", gsea_summary$test), "Pathological", "Clinical")


p1 <- ggplot(gsea_summary, aes(x = NES, y = test)) +
  geom_segment(aes(x = 0, xend = NES, y = test, yend = test, color = `FDR < 0.05`),
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
  labs(x = "NES") +
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

#########################

PCA_summary$test <- factor(PCA_summary$test, levels = rev(gsea_summary$test))

PCA_summary$`FDR < 0.05` <- ifelse(PCA_summary$padj < 0.05, T, F)

PCA_summary$`Disease context` <- ifelse(grepl("CERAD|Braak|AD_dementia|INS_AD|MCI", PCA_summary$test), "AD spectrum", 
                                         ifelse(grepl("Macnair", PCA_summary$test), "Multiple sclerosis",
                                                ifelse(grepl("INS_PSP|INS_PiD", PCA_summary$test), "Primary tauopathy",
                                                       ifelse(grepl("ALS", PCA_summary$test), "ALS-FTD spectrum", 
                                                              ifelse(grepl("PD", PCA_summary$test), "PD-DLB spectrum", "")))))

PCA_summary$`Measure type` <- ifelse(grepl("CERAD|Braak|Rexach|Macnair", PCA_summary$test), "Pathological", "Clinical")


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

svglite(filename = "./analysis/microglia/plots_final/DAM_replication_summary_plot.svg", width = 9.5, height = 10)
plot(p1 + p2)
dev.off()

##################################################
##################################################
##################################################

# plots_final of individual GSEA & PCA

# Green 2024 - PCA/GSEA: Braak
# Macnair 2025 - PCA: all, GSEA: AL
# Rexach 2024 - PCA: all, GSEA: PiD
# Zelic 2025 
# Martirosyan 2024 

DAM_genes <- readRDS("./analysis/microglia/cluster_characterization/final_gene_signatures/DAM_genes_zscore.rds")
pathway <- list(DAM_genes)
names(pathway) <- "DAM"

##################################################
##################################################
##################################################

# Green 2024 - Braak

##################################################

Green_PCA <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Green_2024/PCA_DAM_final.rds")

Green_PCA$PC1 <- Green_PCA$PC1 * -1

p3 <- ggplot(Green_PCA, aes(x = braak_score, y = PC1)) +
  geom_point(aes(colour = braak_score), size = 2) + 
  geom_smooth(method = "lm", colour = "black") +
  scale_color_gradient(low = "#f0cccc", high = "#b30000") +
  theme_bw() +
  labs(x = "Braak score",
       y = "PC1 (17.85%)") + 
  theme(legend.position = "none",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),  
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_x_continuous(breaks = 0:6)

svglite(filename = "./analysis/microglia/plots_final/Green_2024_Braak_PCA.svg", width = 6, height = 3.2)
plot(p3)
dev.off()

##################################################

Green_DE_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Green_2024_braak_deseq.rds")

res_Green_braak <- as.data.frame(results(Green_DE_res, name = "braak_scaled")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Green_braak$t <- sign(res_Green_braak$log2FoldChange) * -log10(res_Green_braak$pvalue)
res_Green_braak <- res_Green_braak %>%
  arrange(desc(t))

ranking_Green_braak <- res_Green_braak$t
names(ranking_Green_braak) <- res_Green_braak$gene

gsea_Green_braak <- fgsea(pathways = pathway, stats = ranking_Green_braak, scoreType = "std")


pd <- plotEnrichmentData(
  pathway = pathway$DAM,,
  stats = ranking_Green_braak)
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
       labs(x="Rank - Braak score", y="Enrichment score"))

svglite(filename = "./analysis/microglia/plots_final/Green_2024_Braak_GSEA.svg", width = 5.8, height = 3.2)
plot(p4)
dev.off()

##################################################
##################################################
##################################################

# Macnair 2025 - AL

##################################################

Macnair_PCA <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Macnair_2025/PCA_DAM_final_WM.rds")

Macnair_PCA$lesion_type <- gsub("WM", "HC-WM", Macnair_PCA$lesion_type)
Macnair_PCA$lesion_type <- gsub("NAHC-WM", "NAWM", Macnair_PCA$lesion_type)

Macnair_PCA$lesion_type <- factor(Macnair_PCA$lesion_type, levels = c("HC-WM", "NAWM", "AL", "CAL", "CIL", "RL"))

p5 <- ggplot(Macnair_PCA, aes(x = lesion_type, y = PC1, fill = lesion_type)) +
  geom_boxplot() + 
  geom_point() + 
  scale_fill_manual(values = disease_colors) +
  theme_bw() +
  labs(x = "Lesion type",
       y = "PC1 (31.25%)") + 
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

svglite(filename = "./analysis/microglia/plots_final/Macnair_2025_PCA.svg", width = 6, height = 3.2)
plot(p5)
dev.off()

##################################################

Macnair_DE_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Macnair_2025_WM_deseq.rds")

res_Macnair_AL <- as.data.frame(results(Macnair_DE_res, name = "lesion_type_AL_vs_WM")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Macnair_AL$t <- sign(res_Macnair_AL$log2FoldChange) * -log10(res_Macnair_AL$pvalue)
res_Macnair_AL <- res_Macnair_AL %>%
  arrange(desc(t))

ranking_Macnair_AL <- res_Macnair_AL$t
names(ranking_Macnair_AL) <- res_Macnair_AL$gene

gsea_Macnair_AL <- fgsea(pathways = pathway, stats = ranking_Macnair_AL, scoreType = "std")


pd <- plotEnrichmentData(
  pathway = pathway$DAM,,
  stats = ranking_Macnair_AL)
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
             labs(x="Rank - AL vs. HC-WM", y="Enrichment score"))

svglite(filename = "./analysis/microglia/plots_final/Macnair_2025_AL_GSEA.svg", width = 5.8, height = 3.2)
plot(p6)
dev.off()

##################################################
##################################################
##################################################

# Rexach 2024 - PiD

##################################################

Rexach_PCA <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Rexach_2024/PCA_DAM_final.rds")

Rexach_PCA$PC1 <- Rexach_PCA$PC1 * -1

p7 <- ggplot(Rexach_PCA, aes(x = npdx1, y = PC1, fill = npdx1)) +
  geom_boxplot() + 
  geom_point() + 
  scale_fill_manual(values = disease_colors) +
  theme_bw() +
  scale_y_continuous(limits = c(-11.5, 12)) + 
  labs(x = "Neuropathological diagnosis",
       y = "PC1 (22.51%)") + 
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

svglite(filename = "./analysis/microglia/plots_final/Rexach_2024_PCA.svg", width = 6, height = 3.2)
plot(p7)
dev.off()

##################################################

Rexach_DE_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Rexach_2024_npdx_deseq.rds")

res_Rexach_PiD <- as.data.frame(results(Rexach_DE_res, name = "npdx1_PiD_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Rexach_PiD$t <- sign(res_Rexach_PiD$log2FoldChange) * -log10(res_Rexach_PiD$pvalue)
res_Rexach_PiD <- res_Rexach_PiD %>%
  arrange(desc(t))

ranking_Rexach_PiD <- res_Rexach_PiD$t
names(ranking_Rexach_PiD) <- res_Rexach_PiD$gene

gsea_Rexach_PiD <- fgsea(pathways = pathway, stats = ranking_Rexach_PiD, scoreType = "std")


pd <- plotEnrichmentData(
  pathway = pathway$DAM,,
  stats = ranking_Rexach_PiD)
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
             labs(x="Rank - PiD vs. HC", y="Enrichment score"))

svglite(filename = "./analysis/microglia/plots_final/Rexach_2024_PiD_GSEA.svg", width = 5.8, height = 3.2)
plot(p8)
dev.off()

##################################################
##################################################
##################################################
# 
# # Zelic 2025 - ALS
# 
# ##################################################
# 
# Zelic_PCA <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Zelic_2025/PCA_DAM_final.rds")
# 
# Zelic_PCA$PC1 <- Zelic_PCA$PC1 * -1
# 
# p9 <- ggplot(Zelic_PCA, aes(x = group, y = PC1, fill = group)) +
#   geom_boxplot() + 
#   geom_point() + 
#   scale_fill_manual(values = disease_colors) +
#   theme_bw() +
#   labs(x = "Diagnosis",
#        y = "PC1 (39.08%)") + 
#   theme(legend.position = "none",
#         text = element_text(family = "Arial"),
#         axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
#         axis.title.y = element_text(face = "bold", size = 14),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12))
# 
# svglite(filename = "./analysis/microglia/plots_final/Zelic_2025_PCA.svg", width = 6, height = 3.2)
# plot(p9)
# dev.off()
# 
# ##################################################
# 
# Zelic_DE_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Zelic_2025_dx_ALS_deseq.rds")
# 
# res_Zelic_ALS <- as.data.frame(results(Zelic_DE_res, name = "group_ALS_vs_HC")) %>%
#   na.omit(.) %>%
#   rownames_to_column(var = "gene")
# 
# res_Zelic_ALS$t <- sign(res_Zelic_ALS$log2FoldChange) * -log10(res_Zelic_ALS$pvalue)
# res_Zelic_ALS <- res_Zelic_ALS %>%
#   arrange(desc(t))
# 
# ranking_Zelic_ALS <- res_Zelic_ALS$t
# names(ranking_Zelic_ALS) <- res_Zelic_ALS$gene
# 
# gsea_Zelic_ALS <- fgsea(pathways = pathway, stats = ranking_Zelic_ALS, scoreType = "std")
# 
# 
# pd <- plotEnrichmentData(
#   pathway = pathway$DAM,,
#   stats = ranking_Zelic_ALS)
# p10 <- with(pd,
#            ggplot(data=curve) +
#              geom_line(aes(x=rank, y=ES), color="green") +
#              geom_segment(data=ticks, mapping=aes(x=rank, y=-spreadES/16, xend=rank, yend=spreadES/16), size=0.2) +
#              geom_hline(yintercept=posES, colour="red", linetype="dashed") +
#              geom_hline(yintercept=negES, colour="red", linetype="dashed") +
#              geom_hline(yintercept=0, colour="black") +
#              theme_bw() +
#              theme(legend.position = "none",
#                    legend.title = element_text(size = 14),
#                    legend.text = element_text(size = 12),  
#                    text = element_text(family = "Arial"),
#                    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
#                    axis.title.y = element_text(face = "bold", size = 14),
#                    axis.text.x = element_text(size = 12),
#                    axis.text.y = element_text(size = 12)) + 
#              labs(x="Rank", y="Enrichment score"))
# 
# svglite(filename = "./analysis/microglia/plots_final/Zelic_2025_ALS_GSEA.svg", width = 5.8, height = 3.2)
# plot(p10)
# dev.off()

##################################################
##################################################
##################################################

# Martirosyan 2024 - PD

##################################################

Martirosyan_PCA <- readRDS("./analysis/microglia/DAM_signature_replication/individual_studies/Martirosyan_2024/PCA_DAM_final.rds")

Martirosyan_PCA$PC1 <- Martirosyan_PCA$PC1 * -1

p11 <- ggplot(Martirosyan_PCA, aes(x = group, y = PC1, fill = group)) +
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

svglite(filename = "./analysis/microglia/plots_final/Martirosyan_2024_PCA.svg", width = 6, height = 3.2)
plot(p11)
dev.off()

##################################################

Martirosyan_DE_res <- readRDS("./analysis/microglia/DAM_signature_replication/DESeq_results/Martirosyan_2024_dx_PD_deseq.rds")

res_Martirosyan_PD <- as.data.frame(results(Martirosyan_DE_res, name = "group_PD_vs_HC")) %>%
  na.omit(.) %>%
  rownames_to_column(var = "gene")

res_Martirosyan_PD$t <- sign(res_Martirosyan_PD$log2FoldChange) * -log10(res_Martirosyan_PD$pvalue)
res_Martirosyan_PD <- res_Martirosyan_PD %>%
  arrange(desc(t))

ranking_Martirosyan_PD <- res_Martirosyan_PD$t
names(ranking_Martirosyan_PD) <- res_Martirosyan_PD$gene

gsea_Martirosyan_PD <- fgsea(pathways = pathway, stats = ranking_Martirosyan_PD, scoreType = "std")


pd <- plotEnrichmentData(
  pathway = pathway$DAM,,
  stats = ranking_Martirosyan_PD)
p12 <- with(pd,
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

svglite(filename = "./analysis/microglia/plots_final/Martirosyan_2024_PD_GSEA.svg", width = 5.8, height = 3.2)
plot(p12)
dev.off()

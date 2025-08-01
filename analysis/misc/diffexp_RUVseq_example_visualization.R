library(tidyverse)
library(ggplot2)
library(svglite)

setwd("/data/ADRD/glia_across_NDDs/")

##################################################
##################################################
##################################################

RUV_ex <- readRDS("./analysis/astrocytes/differential_expression/covar_selection/Mathys_TC_MTG_AD_covar_selection.rds")

##################################################

# elbow plot of residual var

varPart_df <- RUV_ex$varPart_df

p1 <- ggplot(varPart_df, aes(x = RUV_W, y = mean_residual)) + 
  geom_point(aes(size = 3)) +
  geom_line() +
  theme_bw() +
  labs(x = "Component of unwanted variance",
       y = "Mean variance explained\nby the residual",
       title = "Astrocytes, Mathys 2024, MTG, AD vs. HC") +
  theme(text = element_text(family = "Arial"),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

svglite(filename = "./analysis/misc/diffexp_RUVseq_var_explained_by_residual_elbowplot.svg", width = 6.5, height = 3.5)
plot(p1)
dev.off()

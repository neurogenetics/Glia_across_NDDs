library(tidyverse)

folder <- "./Users/acridj/Downloads/Final tables/iMG_DE_results/"#TableS42

# List all CSV files in the folder
csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

# Create a named list to store each data frame
results_list <- lapply(csv_files, read.csv)

# Extract names from filenames and assign as list names
names(results_list) <- tools::file_path_sans_ext(basename(csv_files))

# Combine all data frames into one, adding a 'group' column
combined_df <- do.call(rbind, lapply(names(results_list), function(name) {
   df <- results_list[[name]]
   df$group <- name
   return(df)
}))

# View combined dataframe
head(combined_df)

combined_df<-combined_df%>%
   dplyr::rename(gene=X)

hnDAM <- read.csv("./Users/acridj/Downloads/tableS33.csv")

combined_df%>%
   filter(gene%in%hnDAM$gene)%>%
   mutate(tx_group=case_when(
      group=="Entinostat_0.1μM"~"HDAC",
      group=="Entinostat_1μM"~"HDAC",
      group=="Vorinostat_0.1μM"~"HDAC",
      group=="Vorinostat_1μM"~"HDAC",
      group=="UV-DAMPS_168000"~"cell",
      group=="UV-DAMPS_336000"~"cell",
      group=="UV-DAMPS_672000"~"cell"
   ))%>%
   group_by(tx_group)%>%
   filter(adj.P.Val <0.05 & logFC > 0)%>%
   summarize(n=n_distinct(gene))
#nums: 87/98 hnDAM DE; 12007 of 13479 DE in at least one
#DE in one HDAC|Cell condition: 81H,40c; 45H, 26c up

x<-combined_df%>%
   filter(gene%in%hnDAM$gene)%>%
   mutate(tx_group=case_when(
      group%in%c("Entinostat_0.1μM","Entinostat_1μM","Vorinostat_0.1μM","Vorinostat_1μM")~"HDAC",
      group%in%c("UV-DAMPS_168000","UV-DAMPS_336000","UV-DAMPS_672000")~"cell"
   ))

x<-x%>%
   filter(group%in%c("Entinostat_1μM","Vorinostat_1μM","UV-DAMPS_336000","UV-DAMPS_672000"))#only highest dose

de_hdac<-x%>%filter(adj.P.Val<0.05 & tx_group == "HDAC")

de_hdac_summary <- de_hdac %>%
   group_by(gene) %>%
   summarise(
      n = n(),  # number of entries per gene
      consistent = all(sign(logFC) == sign(logFC[1])),  # TRUE if all same direction
      direction = case_when(
         all(logFC > 0) ~ "up",
         all(logFC < 0) ~ "down",
         TRUE ~ "mixed"
      ),
      .groups = "drop"
   )

de_cell<-x%>%filter(adj.P.Val<0.05 & tx_group == "cell")

de_cell_summary <- de_cell %>%
   group_by(gene) %>%
   summarise(
      n = n(),  # number of entries per gene
      consistent = all(sign(logFC) == sign(logFC[1])),  # TRUE if all same direction
      direction = case_when(
         all(logFC > 0) ~ "up",
         all(logFC < 0) ~ "down",
         TRUE ~ "mixed"
      ),
      .groups = "drop"
   )

####
de_hdac_summary%>%filter(n==2 & consistent==T)%>%group_by(direction)%>%summarise(n=n())

de_cell_summary%>%filter(n==2 & consistent==T)%>%group_by(direction)%>%summarise(n=n())

plot(venn(list("hdac"=de_hdac_summary%>%filter(n==2 & consistent==T)%>%pull(gene),"cell"=de_cell_summary%>%filter(n==2 & consistent==T)%>%pull(gene)))
)


common=intersect(de_hdac_summary%>%filter(n==2 & consistent==T)%>%pull(gene),de_cell_summary%>%filter(n==2 & consistent==T)%>%pull(gene))

test<-combined_df%>%
   filter(gene%in%common)

library(tidyverse)
library(MASS)  # for fitdistr()
library(broom)


x <- read.csv("tableS6.csv")#cell_counts per subcluster (all cells)
y<-x%>%
   group_by(cluster,region)%>%
   summarise(n_cell=sum(n_cells))

# Summary
summary(y$n_cell)
hist(y$n_cell, breaks=30, main="Histogram of Cell Counts", xlab="Cell Counts")

# Fit Negative Binomial
fit_nb <- fitdistr(y$n_cell, "Negative Binomial")
fit_nb

# Plot fitted vs observed
library(ggplot2)
ggplot(y, aes(x = n_cell)) + 
   geom_histogram(binwidth = 1, fill="blue", alpha=0.6) + 
   scale_x_continuous(trans='log1p') +
   labs(title="Histogram of Cell Counts (log scale)", x="log(1 + Cell Counts)")
AIC(zip_model, zinb_model)


###regional anova
y <- y %>%
   mutate(
      major_cell_type = case_when(
         grepl("^Astro", cluster) ~ "astro",
         grepl("^Micro", cluster) ~ "micro",
         grepl("^Oligo", cluster) ~ "oligo",
         TRUE ~ NA_character_
      )
   ) %>%
   filter(!is.na(major_cell_type))



anova_region_effects <- y %>%
   group_by(major_cell_type, cluster) %>%
   do(tidy(aov(n_cell ~ region, data = .))) %>%
   filter(term == "region")  # just keep the region effect
y

region_class_map <- list(
   # cortex = c("PFC", "M1", "TC", "MTG", "OC", "V1"),
   # hippocampus = c("HIP", "HC", "EC"),
   midbrain = c("GPi", "TH", "DMV"),
   # other = c("AnG", "FC")
   cortex = c("PFC", "M1", "TC", "MTG", "OC", "V1","HIP", "HC", "EC","AnG", "FC")
)

region_class_df <- bind_rows(lapply(names(region_class_map), function(class_name) {
   tibble(region = region_class_map[[class_name]], region_class = class_name)
}))

y <- y %>%
   left_join(region_class_df, by = "region")

# anova_regionclass_effects <- y %>%
#    group_by(major_cell_type, cluster) %>%
#    filter(!is.na(region_class)) %>%
#    do(tidy(aov(n_cell ~ region_class, data = .))) %>%
#    filter(term == "region_class")

# Function to run t-test safely
safe_ttest <- function(df) {
   tryCatch({
      tidy(t.test(n_cell ~ region_class, data = df))
   }, error = function(e) {
      tibble(
         estimate = NA,
         estimate1 = NA,
         estimate2 = NA,
         statistic = NA,
         p.value = NA,
         parameter = NA,
         conf.low = NA,
         conf.high = NA,
         method = NA,
         alternative = NA
      )
   })
}

# Run t-test per cluster within each major cell type
ttest_regionclass_effects <- y %>%
   group_by(major_cell_type, cluster) %>%
   filter(!is.na(region_class)) %>%
   filter(n_distinct(region_class) == 2) %>%  # make sure it's binary
   do(safe_ttest(.)) %>%
   ungroup()

anova_region_effects <- anova_region_effects %>%
   group_by(major_cell_type) %>%
   mutate(p_adj = p.adjust(p.value, method = "BH"))


#######
i <- read.csv("tableS11.csv")#cell_counts per subcluster (astro subclust)
library(tidyverse)
k<-i%>%
   group_by(cluster,region)%>%
   summarise(n_cells=sum(n_cells))

k <- k %>%
   mutate(
      # Add the region_class column using case_when
      region_class = case_when(
         region %in% region_class_map$cortex   ~ "cortex",
         region %in% region_class_map$midbrain ~ "midbrain",
         TRUE ~ "other" # This line assigns "other" to any region not in the lists
      )
   ) 
l<-k%>%
   group_by(cluster,region_class)%>%
   summarise(n_cells=sum(n_cells))

summary(l)

write.csv(ttest_regionclass_effects,file = "ttest_regionclass_effects.csv")

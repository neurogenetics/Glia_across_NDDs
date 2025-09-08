library(tidyverse)

setwd("/data/ADRD/glia_across_NDDs")

##################################################

# read in merged DE results files

diffexp_results <- list.files("./analysis/cellculture/iMG_bulk_RNA/tmp_results/", pattern = "\\_df_list.rds$", full.names = T)
diffexp_results <- lapply(diffexp_results, readRDS)
names(diffexp_results)<-c("DA10634","DA10682","DA10725","XR10585")
# diffexp_results <- purrr::flatten(diffexp_results)
flatten_with_names <- function(lst) {
  unlist(
    lapply(names(lst), function(outer_name) {
      inner <- lst[[outer_name]]
      setNames(inner, paste0(outer_name, "_", names(inner)))
    }),
    recursive = FALSE
  )
}

# Use it
diffexp_results<-flatten_with_names(diffexp_results)
names(diffexp_results)

#Now select #XR10585_UV-DAMPS_1.68E+05, XR10585_UV-DAMPS_3.36E+05,XR10585_UV-DAMPS_6.72E+05,XR10585_IFNg_2.50E+00,XR10585_IFNg_1.25E+01,XR10585_IFNg_2.50E+01,DA10634_Entinostat_1.00E-07,DA10634_Entinostat_1.00E-06,DA10634_Vorinostat_1.00E-07,DA10634_Vorinostat_1.00E-06,DA10682_Apilimod_1.00E-08,DA10682_Apilimod_1.00E-07,DA10682_Apilimod_1.00E-06,A10682_Vacuolin_1.00E-07,DA10682_Vacuolin_1.00E-06,A10682_Vacuolin_1.00E-05,DA10682_Ibrutinib_1.00E-08,DA10682_Ibrutinib_1.00E-07,DA10682_Ibrutinib_1.00E-06,DA10682_Spebrutinib_1.00E-07,DA10682_Spebrutinib_1.00E-06,DA10682_Spebrutinib_1.00E-05
selected_names <- c(
  "XR10585_UV-DAMPS_1.68E+05", "XR10585_UV-DAMPS_3.36E+05", "XR10585_UV-DAMPS_6.72E+05",
  "XR10585_IFNg_2.50E+00", "XR10585_IFNg_1.25E+01", "XR10585_IFNg_2.50E+01",
  "DA10634_Entinostat_1.00E-07", "DA10634_Entinostat_1.00E-06",
  "DA10634_Vorinostat_1.00E-07", "DA10634_Vorinostat_1.00E-06",
  "DA10682_Apilimod_1.00E-08", "DA10682_Apilimod_1.00E-07", "DA10682_Apilimod_1.00E-06",
  "DA10682_Vacuolin_1.00E-07", "DA10682_Vacuolin_1.00E-06", "DA10682_Vacuolin_1.00E-05",
  "DA10682_Ibrutinib_1.00E-08", "DA10682_Ibrutinib_1.00E-07", "DA10682_Ibrutinib_1.00E-06",
  "DA10682_Spebrutinib_1.00E-07", "DA10682_Spebrutinib_1.00E-06", "DA10682_Spebrutinib_1.00E-05"
)

long_df <- diffexp_results[selected_names] %>%
  compact() %>%
  imap_dfr(~ .x %>%
             rownames_to_column("gene") %>%
             mutate(comparison = .y))


write.csv(long_df, file = "./analysis/cellculture/iMG_bulk_RNA/differential_expression/cellculture_diffexp_results_for_submission.csv")
saveRDS(diffexp_results[selected_names], file = "./analysis/cellculture/iMG_bulk_RNA/differential_expression/cellculture_diffexp_results_for_submission.rds")

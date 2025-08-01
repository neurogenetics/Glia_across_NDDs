rename_datasets <- function(values) {
  
  # Inner helper function that handles one value
  rename_one <- function(what_to_rename) {
    if (grepl("C9ALS", what_to_rename) & grepl("M1", what_to_rename)) {
      return("Pineda_2024_C9ALS_M1")
    } else if (grepl("ALS", what_to_rename) & grepl("M1", what_to_rename)) {
      return("Pineda_2024_ALS_M1")
    } else if (grepl("C9FTLD", what_to_rename) & grepl("M1", what_to_rename)) {
      return("Pineda_2024_C9FTLD_M1")
    } else if (grepl("FTLD", what_to_rename) & grepl("M1", what_to_rename)) {
      return("Pineda_2024_FTLD_M1")
    } else if (grepl("C9ALS", what_to_rename) & grepl("PFC", what_to_rename)) {
      return("Pineda_2024_C9ALS_PFC")
    } else if (grepl("ALS", what_to_rename) & grepl("PFC", what_to_rename)) {
      return("Pineda_2024_ALS_PFC")
    } else if (grepl("C9FTLD", what_to_rename) & grepl("PFC", what_to_rename)) {
      return("Pineda_2024_C9FTLD_PFC")
    } else if (grepl("FTLD", what_to_rename) & grepl("PFC", what_to_rename)) {
      return("Pineda_2024_FTLD_PFC")
    } else if (grepl("PD", what_to_rename) & grepl("DMV", what_to_rename)) {
      return("NM_2024_PD_DMV")
    } else if (grepl("PD", what_to_rename) & grepl("GPi", what_to_rename)) {
      return("NM_2024_PD_GPi")
    } else if (grepl("PD", what_to_rename) & grepl("M1", what_to_rename)) {
      return("NM_2024_PD_M1")
    } else if (grepl("PD", what_to_rename) & grepl("PFC", what_to_rename)) {
      return("NM_2024_PD_PFC")
    } else if (grepl("PD", what_to_rename) & grepl("V1", what_to_rename)) {
      return("NM_2024_PD_V1")
    } else if (grepl("AD", what_to_rename) & grepl("AnG", what_to_rename)) {
      return("Mathys_2024_AD_AnG")
    } else if (grepl("AD", what_to_rename) & grepl("TH", what_to_rename)) {
      return("Mathys_2024_AD_TH")
    } else if (grepl("AD", what_to_rename) & grepl("EC", what_to_rename)) {
      return("Mathys_2024_AD_EC")
    } else if (grepl("AD", what_to_rename) & grepl("HIP", what_to_rename)) {
      return("Mathys_2024_AD_HIP")
    } else if (grepl("AD", what_to_rename) & grepl("PFC", what_to_rename)) {
      return("Mathys_2024_AD_PFC")
    } else if (grepl("AD", what_to_rename) & grepl("MTG", what_to_rename)) {
      return("Mathys_2024_AD_MTG")
    } else if (grepl("CR", what_to_rename) & grepl("AnG", what_to_rename)) {
      return("Mathys_2024_CR_AnG")
    } else if (grepl("CR", what_to_rename) & grepl("TH", what_to_rename)) {
      return("Mathys_2024_CR_TH")
    } else if (grepl("CR", what_to_rename) & grepl("EC", what_to_rename)) {
      return("Mathys_2024_CR_EC")
    } else if (grepl("CR", what_to_rename) & grepl("HIP", what_to_rename)) {
      return("Mathys_2024_CR_HIP")
    } else if (grepl("CR", what_to_rename) & grepl("PFC", what_to_rename)) {
      return("Mathys_2024_CR_PFC")
    } else if (grepl("CR", what_to_rename) & grepl("MTG", what_to_rename)) {
      return("Mathys_2024_CR_MTG")
    } else if (grepl("GRN", what_to_rename) & grepl("OC", what_to_rename)) {
      return("Gerrits_2022_FTD-GRN_OC")
    } else if (grepl("GRN", what_to_rename) & grepl("TC", what_to_rename)) {
      return("Gerrits_2022_FTD-GRN_TC")
    } else if (grepl("GRN", what_to_rename) & grepl("FC", what_to_rename)) {
      return("Gerrits_2022_FTD-GRN_FC")
    } else {
      return(what_to_rename)  # return original if no match
    }
  }
  
  # Apply the helper function to each item
  vapply(values, rename_one, character(1))
}

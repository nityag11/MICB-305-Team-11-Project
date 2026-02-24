# Load in required packages + metadata
library(dplyr)
library(tidyverse)

ms_metadata <- read.delim("ms_metadata_corrected.tsv")

# Check for unique values
unique(ms_metadata$administration)

# Group Subcutaneous and Infused into Non-Oral
ms_treatments_grouped <- ms_metadata %>%
  mutate(treatments_grouped = case_when(
    administration == "Subcutaneous" ~ "Non-Oral",
    administration == "Infused" ~ 'Non-Oral',
    administration == "Untreated" ~ "Untreated",
    administration == "Control" ~ "Control",
    administration == "Oral" ~ "Oral"
  ))
  
ms_treatments_grouped <- ms_treatments_grouped %>%
  relocate(treatments_grouped, .after = administration)

#Convert probiotics to categorical Yes/No
ms_curated <- ms_treatments_grouped %>%
  mutate(probiotics = case_when(
    probiotics == 0 ~ "No",
    probiotics == 1 ~ 'Yes'
  ))

#Save as .tsv
readr::write_tsv(ms_curated, "ms_metadata_curated.tsv")

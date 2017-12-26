
library(tidyverse)
np = readRDS("scratch/processed/cleaned_np_civic_data.RDS") %>%
  filter(sample_id == "A")

cps = readRDS("scratch/processed/cps_civic_full_edited.RDS") %>%
  sample_n(3000)
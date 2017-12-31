
library(tidyverse)
np = readRDS("scratch/processed/cleaned_np_civic_data.RDS") %>%
  filter(sample_id == "A") %>%
  sample_n(50)

cps = readRDS("scratch/processed/cps_civic_full_edited.RDS") %>%
  sample_n(52)

weights = fpbb_synth_pops(cps$pwsrwgt)

x_vars = c("age", "sex", "racethn")
y_vars = c("y_always_vote_local", "y_civic_assoc_yes")

bart_params = list(ndpost = 25, 
                   keeptrainfits = TRUE,
                   mc.cores = 1, 
                   ntree = 50)


bestimate(np, cps, y_vars, x_vars, weights)



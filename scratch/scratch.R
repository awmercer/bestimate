
library(tidyverse)
np = readRDS("scratch/processed/cleaned_np_civic_data.RDS") %>%
  filter(sample_id == "A") %>%
  sample_n(100)

cps = readRDS("scratch/processed/cps_civic_full_edited.RDS") %>%
  sample_n(200)

weights = fpbb_synth_pops(cps$pwsrwgt, L = 1, N = nrow(cps) * 100)

x_vars = c("age", "sex", "racethn", "educcat", "fcregion")
y_vars = c("y_always_vote_local")

bart_params = list(ndpost = 25, 
                   keeptrainfits = TRUE,
                   mc.cores = 3, 
                   ntree = 50)


x = bestimate(np, cps, y_vars, x_vars, weights, pred_subsample_size = 150, model_params = bart_params)



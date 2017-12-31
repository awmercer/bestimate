
library(tidyverse)
np = readRDS("scratch/processed/cleaned_np_civic_data.RDS") %>%
  filter(sample_id == "A") %>%
  sample_n(100)

cps = readRDS("scratch/processed/cps_civic_full_edited.RDS") %>%
  sample_n(200)

sp_weights = fpbb_synth_pops(cps$pwsrwgt, L = 1, N = nrow(cps) * 100)

x_vars = c("age", "sex", "racethn", "educcat", "fcregion")
y_vars = str_subset(names(np), "^y_")

bart_params = list(ndpost = 25, 
                   keeptrainfits = TRUE,
                   mc.cores = 3, 
                   ntree = 50)


x = bestimate(samp = np, 
              ref = cps, 
              y_var_names = y_vars,
              x_var_names = x_vars,
              sp_wts = sp_weights, 
              propensity = TRUE, 
              prediction = TRUE, 
              double_robust = TRUE, 
              dr_propensity_transform = odds, 
              decompose_bias = TRUE, 
              pred_subsample_size = 100, 
              posterior_draws = 1000, 
              model_params = list(keeptrainfits = TRUE, 
                                  mc.cores = 1, 
                                  ntree = 50),
              seed = 12345)

              
              
              



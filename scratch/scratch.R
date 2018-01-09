
library(tidyverse)
library(stringr)
np = readRDS("~/fpbb-inference/data/cleaned/cleaned_np_civic_data.RDS") %>%
  filter(sample_id == "A") %>%
  sample_n(100)

cps = readRDS("~/fpbb-inference/data/cleaned/cps_civic_full_edited.RDS") %>%
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
              propensity = FALSE, 
              prediction = TRUE, 
              double_robust_reg = FALSE, 
              double_robust_wt = FALSE,
              dr_propensity_transform = odds, 
              decompose_bias = FALSE, 
              quality_measures = FALSE,
              pred_subsample_size = 100, 
              posterior_draws = 100, 
              model_params = list(keeptrainfits = TRUE, 
                                  mc.cores = 1, 
                                  ntree = 50),
              seed = 12345)

              
              
              



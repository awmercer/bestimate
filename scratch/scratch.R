
library(tidyverse)
library(stringr)
np = readRDS("~/fpbb-inference/data/cleaned/cleaned_np_civic_data.RDS") %>%
  filter(sample_id == "A")

cps = readRDS("~/fpbb-inference/data/cleaned/cps_civic_full_edited.RDS") %>%
  sample_n(5000)

sp_weights = fpbb_synth_pops(cps$pwsrwgt, L = 1, N = nrow(cps) * 100)

x_vars = c("age", "sex", "racethn", "educcat", "fcregion")
y_vars = str_subset(names(np), "^y_")

bart_params = list(ndpost = 25, 
                   keeptrainfits = TRUE,
                   mc.cores = 3, 
                   ntree = 50)

x_samp = np[, x_vars]
x_ref = cps[, x_vars]

p = imap(sp_weights, ~estimate_propensities(sp_wt = .x, sp_id = .y, x_samp = x_samp, x_ref = x_ref, num_replicates = 2,
                                            bart_params = bart_params))

ref = cps
samp = np

comb = bind_rows(ref %>% mutate(smp = 0), 
                 samp %>% mutate(smp = 1)) %>%
  select(x_vars, smp)


# p = balanced_bagged_bart(x.train = comb[, x_vars], 
#                          y.train = comb[["smp"]], 
#                          num_fits = 2, mc.cores = 2)
# 
# x = predict(p, newdata = samp)

x = bestimate(samp = np, 
              ref = cps, 
              y_var_names = y_vars,
              x_var_names = x_vars,
              sp_wts = sp_weights, 
              propensity = T, 
              prediction = T, 
              double_robust_reg = T, 
              double_robust_wt = T,
              dr_propensity_transform = log, 
              quality_measures = T,
              pred_subsample_size = 100, 
              posterior_draws = 100, 
              bart_params = list(keeptrainfits = TRUE, 
                                  mc.cores = 3, 
                                  ntree = 50))

              
              
              




library(tidyverse)
library(stringr)
np = readRDS("~/fpbb-inference/data/cleaned/cleaned_np_civic_data.RDS") %>%
  filter(sample_id == "A")

cps = readRDS("~/fpbb-inference/data/cleaned/cps_civic_full_edited.RDS") 

sp_weights = fpbb_synth_pops(cps$pwsrwgt, L = 1, N = nrow(cps) * 100)

x_vars = c("age", "sex", "racethn", "educcat", "fcregion")
y_vars = str_subset(names(np), "^y_")


b = pbart2(np[, x_vars], np$y_talk_neighbor_weekly)

b = balanced_bagged_bart(x.train = np[, x_vars], y.train = np$y_always_vote_local, num_fits = 1)

reps1 = bestimate(samp = np, 
                   ref = cps, 
                   y_var_names = y_vars[1],
                   x_var_names = x_vars,
                   sp_wts = sp_weights, 
                   propensity = T, 
                   prediction = F, 
                   double_robust_reg = F, 
                   double_robust_wt = F,
                   dr_propensity_transform = log, 
                   propensity_replicates = 1,
                   quality_measures = F,
                   posterior_draws = 1000, 
                   bart_params = list(mc.cores = 10, 
                                      ntree = 50))
reps10 = bestimate(samp = np, 
              ref = cps, 
              y_var_names = y_vars[1],
              x_var_names = x_vars,
              sp_wts = sp_weights, 
              propensity = T, 
              prediction = F, 
              double_robust_reg = F, 
              double_robust_wt = F,
              dr_propensity_transform = log, 
              propensity_replicates = 10,
              quality_measures = F,
              posterior_draws = 1000, 
              bart_params = list(mc.cores = 10, 
                                  ntree = 50))

reps20 = bestimate(samp = np, 
                   ref = cps, 
                   y_var_names = y_vars[1],
                   x_var_names = x_vars,
                   sp_wts = sp_weights, 
                   propensity = T, 
                   prediction = F, 
                   double_robust_reg = F, 
                   double_robust_wt = F,
                   dr_propensity_transform = log, 
                   propensity_replicates = 20,
                   quality_measures = F,
                   posterior_draws = 1000, 
                   bart_params = list(mc.cores = 10, 
                                      ntree = 50))              
              

var(reps1$y_bar_propwt)

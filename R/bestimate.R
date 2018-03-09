

#' Bayesian finite population with BART and FPBB
#'
#' @param samp \code{data.frame} containing the survey sample data.
#' @param ref  \code{data.frame} containing the reference sample.
#' @param y_var_names Character vector containing names of the outcome
#'     variables.
#' @param x_var_names Character vector containing names of the covariates
#' @param sp_wt_frame \code{data.frame} containing one column of weights
#'     for each synthetic population.
#' @param posterior_draws Number of posterior draws
#'
#' @return Good question
#' @export
bestimate = function(samp,
                     ref,
                     y_var_names,
                     x_var_names,
                     sp_wts = NULL,
                     propensity = TRUE,
                     prediction = TRUE,
                     double_robust_wt = TRUE,
                     double_robust_reg = TRUE,
                     dr_propensity_transform = log,
                     propensity_replicates = 1,
                     quality_measures = TRUE,
                     posterior_draws = 1000,
                     bart_params = list(mc.cores = 1,
                                         ntree = 50)) {
  
  t_start = proc.time()
  
  # If no sp_wts supplied set to vector of 1's. Assumes that the full
  # syntehtic population was passed in as ref
  if (is.null(sp_wts)) {
    sp_wts = tibble(pop = rep(1, nrow(ref)))
  }
  
  # Override required BART parameters
  bart_params$keeptrainfits = TRUE
  bart_params$ndpost = posterior_draws
  
  # Calculate actual number of posterior draws that BART will
  # perform if more than 1 core is used. This will be the
  # lowest multiple of mc.cores that is >= ndpost
  if (!is.null(bart_params$mc.cores)) {
    posterior_draws = ceiling(posterior_draws/bart_params$mc.cores) * bart_params$mc.cores
  }
  
  tr_name = deparse(substitute(dr_propensity_transform))
  
  y_var_names = set_names(y_var_names)
  n_samp = nrow(samp)
  n_ref = nrow(ref)
  
  y_samp = select(samp, y_var_names)
  y_ref = select(ref, y_var_names)
  x_samp = select(samp, x_var_names)
  x_ref = select(ref, x_var_names)
  
  sp_names = names(sp_wts) %>% set_names()
  
  
  # Fit propensity model for each synthetic population
  if (propensity | double_robust_reg | double_robust_wt |  quality_measures) {
    cat("Fitting propensity models:\n")

    # Fit a propensity model to each subsample
    propensities = sp_wts %>%
      imap(
        ~ estimate_propensities(
          sp_wt = .x,
          sp_id = .y,
          x_samp = x_samp,
          x_ref = x_ref,
          bart_params = bart_params
        )
      ) %>%
      transpose()
    
    # Create weights for common support and no common support
    cs_wts = map2(propensities$ref_phi_posterior, sp_wts,
                  function(phi, wt) {
                    map_dfc(phi, ~ . * wt)
                  })
    
    ncs_wts = map2(propensities$ref_phi_posterior, sp_wts,
                   function(phi, wt) {
                     map_dfc(abs(phi - 1), ~ . * wt)
                   })
    
    ref_wts = list(tot = sp_wts, cs = cs_wts, ncs = ncs_wts)
    
    
    # Bayesian bootstrap weights for unweighted sample data
    bb_wts = matrix(rexp(n_samp * posterior_draws, 1),
                    nrow = n_samp,
                    byrow = FALSE)
    bb_wts = bb_wts / colSums(bb_wts)
    
  }
  

  # Generate all estimates for each y_var_name
  estimate_posteriors = y_var_names %>%
    map(function(y_var_name) {
      
      t1 = proc.time()
      cat(sprintf("Starting: %s\n", y_var_name))
      
      res = list()
      y_obs_samp = y_samp[[y_var_name]]
      
      # Primary estimates: propensity, prediction, double-robust
      if (propensity) {
        cat("  Propensity weighted estimates\n")
        res$y_bar_propwt = map(
          propensities$samp_propensity_wts,
          ~ weighted_mean_posterior(y = y_obs_samp, 
                                    w = .x)
        )
      }
      
      if (prediction | double_robust_wt ) {
        cat("  Fitting prediction models\n")
        pred_fit = run_bart_fit(
          y_train = y_obs_samp,
          x_train = x_samp,
          x_test = x_ref,
          bart_params = bart_params
        ) %>% (function(l) {
          list(samp_posterior = l$train_posterior,
               ref_posterior = l$test_posterior)
        })
      }
      if (prediction | double_robust_wt) {
          res$y_bar_pred = map(sp_wts,
                             ~ weighted_mean_posterior(y = pred_fit$ref_posterior,
                                                       w = .x))
      }
      
      if (double_robust_wt) {
        
        term1 =  res$y_bar_pred
        
        residual = map_dfc(pred_fit$samp_posterior, 
                        ~ y_obs_samp - .x
                        )
        term2 = map(propensities$samp_propensity_wts,
                    ~ weighted_mean_posterior(y = residual, 
                                              w = .x))
       
        res$y_bar_drrbc = map2(term1, term2, ~ .x + .y)
        
        if (quality_measures) {
          res$mean_drrbc_term2 = term2 
        }
      }
      
      if (double_robust_reg) {
        cat("  Fitting double-robust models\n")
        tr_prop_name = sprintf("propensity_%s", tr_name)
        
        samp_tr_prop = map(
          propensities$samp_mean_propensity,
          ~ tibble(!!tr_prop_name := dr_propensity_transform(.x))
        )
        
        ref_tr_prop = map(
          propensities$ref_mean_propensity,
          ~ tibble(!!tr_prop_name := dr_propensity_transform(.x))
        )
        
        dr_fits = pmap(
          .f = fit_double_robust,
          .l = list(
            propensity_train = samp_tr_prop,
            propensity_test = ref_tr_prop
          ),
          y_train = y_obs_samp,
          x_train = x_samp,
          x_test = x_ref,
          bart_params = bart_params
        ) %>% map(function(l) {
          list(samp_posterior = l$train_posterior,
               ref_posterior = l$test_posterior)
        }) %>%
          transpose()
        
        res$y_bar_drpsc = map2(dr_fits$ref_posterior,
                            sp_wts,
                            ~ weighted_mean_posterior(y = .x, w = .y))
      }
      

      y_obs_ref = y_ref[[y_var_name]]
      

      
      if (quality_measures) {
        # y_bar_obs_samp (Bayesian bootstrap)
        y_bar_samp_bayesboot = weighted_mean_posterior(y = y_obs_samp,
                                                       w = bb_wts)
        res$y_bar_samp_bayesboot = map(sp_names, ~ y_bar_samp_bayesboot)
        
        
        # Note, because y_obs_ref is a single vector, we
        # duplicate it posterior_draws times so that it works
        # with ref_mean_components. 
        tmp = rep(as.tibble(y_obs_ref), posterior_draws) %>%
          bind_cols()
        res$y_bar_obs_ref = ref_mean_components(tmp, ref_wts, name = "y_bar_obs_ref")
        rm(tmp)
        
        # pct_common_support
        res$pct_common_support = map2(cs_wts, sp_wts, ~ colSums(.x) / sum(.y))
        
        # min_propensity
        res$min_propensity = propensities$samp_propensity_mins

        if (prediction | double_robust_wt) {
          ## BART Model error on reference sample
          pred_model_errors = map_dfc(pred_fit$ref_posterior,
                                      ~ .x - y_obs_ref)
          
          res$pred_model_bias = ref_mean_components(pred_model_errors,
                                                    ref_wts,
                                                    "pred_model_bias")
          
          res$pred_model_rmse = pred_model_errors %>%
            map_dfc( ~ .x ^ 2) %>%
            ref_mean_components(ref_wts, "pred_model_rmse") %>%
            map(sqrt)
        }
        
        if (double_robust_reg) {
          dr_model_errors = map(dr_fits$ref_posterior,
                                ~ map_dfc(.x, ~ .x - y_obs_ref))
          
          res$dr_model_bias = ref_mean_components(dr_model_errors, ref_wts, "dr_model_bias")
          
          res$dr_model_rmse = dr_model_errors %>%
            map( ~ .x ^ 2) %>%
            ref_mean_components(ref_wts, "dr_model_rmse") %>%
            map(sqrt)
        }

      }
      cat("Combining results\n")
      
      out = imap(res, combine_results) %>%
        bind_cols()
      t2 = proc.time() - t1
      cat(sprintf("Finished %s: -- %.1f seconds\n", y_var_name, t2[[3]]))
      
      # Append sp and draw ids
      map(sp_names, ~tibble(sp = .x, draw = seq_len(posterior_draws))) %>%
        bind_rows() %>%
        bind_cols(out)
      
    }) %>% # End of anonymous function
    bind_rows(.id = "y_var")
  
  t_end = proc.time() - t_start
  cat(sprintf("Finished everything: %.1f seconds\n", t_end[[3]]))
  
  return(estimate_posteriors)
  
}


combine_results = function(l, name) {
  if (is.data.frame(l[[1]])) {
    return(bind_rows(l))
  } else {
    tibble(!!name := unlist(l))
  }
}

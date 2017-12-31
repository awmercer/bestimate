

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
                     sp_wts,
                     propensity = TRUE,
                     prediction = TRUE,
                     double_robust = TRUE,
                     dr_propensity_transform = odds,
                     decompose_bias = TRUE,
                     pred_subsample_size = nrow(ref),
                     posterior_draws = 25,
                     model_params = list(keeptrainfits = TRUE,
                                         mc.cores = 1,
                                         ntree = 50),
                     seed = 12345) {
  model_params$ndpost = posterior_draws
  
  tr_name = deparse(substitute(dr_propensity_transform))
  
  y_var_names = set_names(y_var_names)
  n_samp = nrow(samp)
  
  y_samp = select(samp, y_var_names)
  y_ref = select(ref, y_var_names)
  x_samp = select(samp, x_var_names)
  x_ref = select(ref, x_var_names)
  
  sp_names = names(sp_wts) %>% set_names()
  
  
  # Fit propensity model for each synthetic population
  if (propensity | double_robust | decompose_bias) {
    cat("Fitting propensity models:\n")
    
    
    # Create balanced subsamples for estimating
    # propensities
    ref_subsamples_prop = sp_subsample(
      synth_pops = sp_wts,
      size = n_samp,
      as_weights = TRUE,
      seed = seed
    )
    
    # Fit a propensity model to each subsample
    propensities = ref_subsamples_prop %>%
      imap(
        ~ estimate_propensities(
          subsamp_ids = .x,
          sp_id = .y,
          x_samp = x_samp,
          x_ref = x_ref,
          model_params = model_params,
          seed = seed
        )
      ) %>%
      transpose()
  }
  
  if (decompose_bias) {
    # Create random subsamples for unconfounded bias decomposition models
    if (!is.null(seed))
      set.seed(seed)
    
    ref_subsamples_pred = sp_subsample(
      synth_pops = sp_wts,
      size = pred_subsample_size,
      as_weights = TRUE,
      seed = seed
    )
    
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
    if (!is.null(seed))
      set.seed(seed)
    bb_wts = matrix(rexp(n_samp * posterior_draws, 1),
                    nrow = n_samp,
                    byrow = FALSE)
    bb_wts = bb_wts / colSums(bb_wts)
    
  }
  
  
  # Generate all estimates for each y_var_name
  y_var_names %>%
    map(function(y_var_name) {
      res = list()
      y_obs_samp = y_samp[[y_var_name]]
      
      # Primary estimates: propensity, prediction, double-robust
      if (propensity) {
        res$y_bar_propwt = map(
          propensities$samp_propensity_wts,
          ~ weighted_mean_posterior(y = y_obs_samp, 
                                    w = .x)
        )
      }
      
      if (prediction | decompose_bias) {
        pred_fit = run_bart_fit(
          y_train = y_obs_samp,
          x_train = x_samp,
          x_test = x_ref,
          bart_params = model_params,
          seed = seed
        )
        res$y_bar_pred = map(sp_wts,
                             ~ weighted_mean_posterior(y = pred_fit$test_posterior,
                                                       w = .x))
      }
      
      if (double_robust) {
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
          model_params = model_params,
          seed = seed
        ) %>%
          transpose()
        
        res$y_bar_dr = map2(dr_fits$test_posterior,
                            sp_wts,
                            ~ weighted_mean_posterior(y = .x, w = .y))
      }
      
      # Secondary estimates for bias decomposition
      if (decompose_bias) {
        y_obs_ref = y_ref[[y_var_name]]
        
        pred_fits_unconfounded = map(
          ref_subsamples_pred,
          ~ run_bart_fit(
            y_train = y_obs_ref[.x],
            x_train = slice(x_ref, .x),
            x_test = x_samp,
            bart_params = model_params,
            seed = seed
          )
        ) %>% transpose()
        
        if (double_robust) {
          dr_fits_unconfounded = pmap(list(ref_subsamples_pred,
                                           ref_tr_prop,
                                           samp_tr_prop),
                                      function(s, pr_train, pr_test) {
                                        fit_double_robust(
                                          y_train = y_obs_ref[s],
                                          x_train = slice(x_ref, s),
                                          propensity_train = slice(pr_train, s),
                                          x_test = x_samp,
                                          propensity_test = pr_test,
                                          model_params = model_params,
                                          seed = seed
                                        )
                                      }) %>%
            transpose()
        }
        
        
        # To Create - conditional means
        # y_bar_obs_samp (Bayesian bootstrap)
        y_bar_samp_bayesboot = weighted_mean_posterior(y = y_obs_samp,
                                                       w = bb_wts)
        res$y_bar_samp_bayesboot = map(sp_names, ~ y_bar_samp_bayesboot)
        
        res$y_bar_samp_confounded = list(colMeans(pred_fit$train_posterior)) %>%
          map2(sp_names, ~ .x) %>%
          set_names(sp_names)
        
        res$y_bar_samp_unconfounded = map(pred_fits_unconfounded$test_posterior,
                                          colMeans)
        
        
        # Note, because y_obs_ref is a single vector, we
        # duplicate it posterior_draws times so that it works
        # with ref_mean_components
        tmp = rep(as.tibble(y_obs_ref), posterior_draws) %>%
          bind_cols()
        res$y_bar_obs_ref = ref_mean_components(tmp, ref_wts, name = "y_bar_obs_ref")
        rm(tmp)
        
        if (propensity) {
          # y_bar_propwt_confounded
          res$y_bar_propwt_confounded = map(
            propensities$samp_propensity_wts,
            ~ weighted_mean_posterior(y = pred_fit$train_posterior,
                                      w = .x)
          )
          
          # y_bar_propwt_unconfounded
          res$y_bar_propwt_unconfounded = map2(
            pred_fits_unconfounded$test_posterior,
            propensities$samp_propensity_wts,
            weighted_mean_posterior
          )
        }
        
        if (prediction) {
          # y_bar_pred_confounded (tot, cs, ncs)
          res$y_bar_pred_confounded = ref_mean_components(
            y_pos = pred_fit$test_posterior,
            wts = ref_wts,
            name = "y_bar_pred_confounded"
          )
          
          # y_bar_pred_unconfounded (tot, cs, ncs)
          res$y_bar_pred_unconfounded = ref_mean_components(
            y_pos = pred_fits_unconfounded$train_posterior,
            wts = ref_wts,
            name = "y_bar_pred_unconfounded"
          )
          
        }
        if (double_robust) {
          # y_bar_dr_confounded (tot, cs, ncs)
          res$y_bar_dr_confounded = ref_mean_components(
            y_pos = dr_fits$test_posterior,
            wts = ref_wts,
            name = "y_bar_dr_confounded"
          )
          
          # y_bar_dr_unconfounded (tot, cs, ncs)
          res$y_bar_dr_unconfounded = ref_mean_components(
            y_pos = dr_fits_unconfounded$train_posterior,
            wts = ref_wts,
            name = "y_bar_dr_unconfounded"
          )
        }
        # To Create - other quality measures
        
        delta_i_samp_pred = map(pred_fits_unconfounded$test_posterior,
                                ~ pred_fit$train_posterior - .x)
        
        res$delta_exch_unwt = map(delta_i_samp_pred, colMeans)
        res$delta_exch_propwt = map2(
          delta_i_samp_pred,
          propensities$samp_propensity_wts,
          weighted_mean_posterior
        )
        
        rm(delta_i_samp_pred)
        
        if (prediction) {
          res$delta_exch_pred = map(pred_fits_unconfounded$train_posterior,
                                    ~ pred_fit$test_posterior  - .x) %>%
            ref_mean_components(ref_wts, "delta_exch_pred")
        }
        
        if (double_robust) {
          res$delta_exch_dr = map2(dr_fits$test_posterior,
                                   dr_fits_unconfounded$train_posterior,
                                   ~ .x - .y) %>%
            ref_mean_components(ref_wts, "delta_exch_dr")
        }
        
        # pct_common_support
        res$pct_common_support = map2(cs_wts, sp_wts, ~ colSums(.x) / sum(.y))
        
        # min_propensity
        res$min_propensity = propensities$samp_propensity_mins
        
        ## BART Model error on reference sample
        pred_model_errors = map_dfc(pred_fit$test_posterior,
                                    ~ .x - y_obs_ref)
        
        if (prediction) {
          res$pred_model_bias = ref_mean_components(pred_model_errors,
                                                    ref_wts,
                                                    "pred_model_bias")
          
          res$pred_model_rmse = pred_model_errors %>%
            map_dfc( ~ .x ^ 2) %>%
            ref_mean_components(ref_wts, "pred_model_rmse") %>%
            map(sqrt)
        }
        
        if (double_robust) {
          dr_model_errors = map(dr_fits$test_posterior,
                                ~ map_dfc(.x, ~ .x - y_obs_ref))
          
          res$dr_model_bias = ref_mean_components(dr_model_errors, ref_wts, "dr_model_bias")
          
          res$dr_model_rmse = dr_model_errors %>%
            map( ~ .x ^ 2) %>%
            ref_mean_components(ref_wts, "dr_model_rmse") %>%
            map(sqrt)
        }
        
        
        
      }
      imap(res, combine_results) %>%
        bind_cols()
    }) %>% # End of anonymous function
    bind_rows(.id = "y_var")
  
}


combine_results = function(l, name) {
  if (is.data.frame(l[[1]])) {
    return(bind_rows(l))
  } else {
    tibble(!!name := unlist(l))
  }
}

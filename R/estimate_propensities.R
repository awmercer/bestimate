

#' Title
#'
#' @param subsamp_ids 
#' @param sp_id 
#' @param x_samp 
#' @param x_ref 
#' @param model_params 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
estimate_propensities = function(subsamp_ids, 
                            sp_id, 
                            x_samp, 
                            x_ref, 
                            model_params,
                            seed
                            ) {
  t = proc.time()
  cat(sprintf("    %s ", sp_id))
  
  n_samp = nrow(x_samp)
  # Combine sample and subsample from reference data
  comb = bind_rows(x_samp, 
                   slice(x_ref, subsamp_ids))
  
  origin = c(rep(1, n_samp), rep(0, n_samp))
  
  # Override keeptrainfits and set to TRUE if set to FALSE
  model_params$keeptrainfits = TRUE
  prop_posteriors = run_bart_fit(x_train = comb,
                          y_train = origin, 
                          x_test = x_ref, 
                          model_params,
                          seed=seed)
  
  samp_propensity_posterior = prop_posteriors$train_posterior[1:n_samp, ]
  ref_propensity_posterior = prop_posteriors$test_posterior

  samp_min = apply(samp_propensity_posterior, 2, min)
  
  ref_phi_posterior = seq_along(samp_min) %>%
    set_names(sprintf("draw_%s", .)) %>%
    map(~as.numeric(ref_propensity_posterior[, .] >= samp_min[.])) %>%
    bind_cols()

  t2 = proc.time() - t
  cat(sprintf("%.1f\n", t2[[3]]))
  
  list(ref_mean_propensity = rowMeans(ref_propensity_posterior),
       samp_mean_propensity = rowMeans(samp_propensity_posterior),
       samp_propensity_posterior = samp_propensity_posterior,
       samp_propensity_wts = prop2weight(samp_propensity_posterior, 
                                         p_denom = TRUE),
       ref_phi_posterior = ref_phi_posterior,
       samp_propensity_mins = samp_min,
       propensity_fit = prop_posteriors$fit)
}







#' Title
#'
#' @param subsamp_ids 
#' @param sp_id 
#' @param x_samp 
#' @param x_ref 
#' @param bart_params 
#'
#' @return
#' @export
#'
#' @examples
estimate_propensities = function(
                            sp_wt,
                            sp_id, 
                            x_samp, 
                            x_ref, 
                            bart_params,
                            num_replicates=10
                            ) {
  t = proc.time()
  cat(sprintf("    %s ", sp_id))
  
  n_samp = nrow(x_samp)
  n_ref = nrow(x_ref)
  # Combine sample and subsample from reference data
  comb = bind_rows(x_samp, 
                   x_ref)
  
  origin = c(rep(1, n_samp), rep(0, n_ref))
  comb_wts = c(rep(1, n_samp), sp_wt)
  
  prop_fits = balanced_bagged_bart(x.train = comb, 
                                   y.train = origin, 
                                   num_fits = num_replicates, 
                                   weights = comb_wts, 
                                   bart_params = bart_params)
  
  mc.cores = ifelse(is.null(bart_params$mc.cores), 1, bart_params$mc.cores)

  samp_propensity_posterior = predict(prop_fits, newdata = x_samp, mc.cores = mc.cores)
  ref_propensity_posterior = predict(prop_fits, newdata = x_ref, mc.cores = mc.cores)

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
       propensity_fit = prop_fits)
}





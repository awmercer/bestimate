
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
#'
bestimate = function(samp,
                     ref,
                     y_var_names,
                     x_var_names,
                     sp_wts,
                     method = c("propensity", "prediction", "dr"),
                     posterior_draws = 100,
                     mc.cores = 1,
                     ... # Additional options passed to BART
                     ) {
  
  # Check to make sure that the X vars exist and are compatible across
  # samp and ref
  # Check to make sure y_vars exist 
  
  x_samp = select(samp, x_var_names) 
  x_ref = select(ref, x_var_names) 
  
  res = list(predictions = NULL,
             propensities = NULL,
             dr_predictions = NULL)

  # Fit prediction model for each y_variable specified
  if ("prediction" %in% method) {
    res[["predictions"]] = y_var_names %>%
      set_names(.) %>%
      map( function(y_var, ...) {
        pbart2(x.train = x_samp, 
               y.train = samp[[y_var]],
               x.test = x_ref,
               ndpost = posterior_draws, 
               mc.cores = mc.cores,
               ...
                  )}, ...) %>%
      map("yhat.test") %>% # Extract the posterior for the ref sample
      map(pnorm)        # Convert to probabilities
  }
  
  # Fit propensity model for each synthetic population
  if ("propensity" %in% method | "dr" %in% method) {
    
    # Get a random subsample from each synthetic population
    # of a size equal to that of the survey sample
    n_samp = nrow(samp)
    
    ref_subsamples = sp_wts %>%
      map(expand_sp_weights) %>%
      map(sample, size=n_samp, replace=FALSE) 
    

    # Fit a propensity model to each subsample
    res[["propensities"]] = ref_subsamples %>%
      map(function(subsamp_ids, ...) {
        
        # Combine sample and subsample from reference data
        comb = bind_rows(x_samp, 
                         slice(x_ref, subsamp_ids))
        
        # Fit propensity model to balanced dataset
        prop_fit = pbart2(x.train = comb, 
                          y.train = c(rep(1, n_samp), rep(0, n_samp)),
                          x.test = x_ref, 
                          ndpost = posterior_draws, 
                          keeptrainfits = TRUE,
                          mc.cores = mc.cores,
                          ...)
        
        samp_p = prop_fit$yhat.train.mean[1:n_samp] %>% pnorm()
        ref_p = prop_fit$yhat.test.mean %>% pnorm()
        
        # Return only the mean propensities for ref and sample plus
        # The full posterior distribution for the sample
        list(ref_mean_propensity = ref_p,
             samp_mean_propensity = samp_p,
             samp_propensity_posterior = prop_fit$yhat.train[, 1:n_samp] %>% 
               pnorm() %>% t())
      }, ...)
  }
  
  if ("dr" %in% method) {
    # For each synthetic population, fot a prediction model for each Y
    # The prediction model is fit using the odds of the mean propensity
  res[["dr_predictions"]] = cross(list(y_var = y_var_names, 
                                       sp = names(sp_wts))) %>%
    transpose() %>%
    pmap(function(y_var, sp, ...) {
      p = res[["propensities"]][[sp]]
      
      x_ref_p = bind_cols(x_ref, prop__ = odds(p$ref_mean_propensity))
      
      x_samp_p = bind_cols(x_samp, prop__ = odds(p$samp_mean_propensity))
      
      dr_fit = pbart2(x.train = x_samp_p, 
                      y.train = samp[[y_var]], 
                      x.test = x_ref_p, 
                      ndpost = posterior_draws,
                      mc.cores = mc.cores,
                      ...
                      )
      
      
    }, ...) %>%
    map("yhat.test") %>%
    map(pnorm)
  }
  
  
  # Fit propensity model to balanced samples
  # If prediction:
  # Fit a prediction model for each outcome variable
  # If DR
  # Append the mean propensity weight to the reference sample
  # Fit a DR model for each outcome variable/synthetic pop
  
  return(res)
}


### ----------------------------------------------------------



# Utility function to multiply a matrix of 
# y values by a matrix of weights and return
# the weighted posterior distribution for y
#' Calculate the posterior distribution for Y over a matrix of weights
#'
#' @param y A data frame or matrix containing one or more columns of Y variables
#' @param w A data frame or matrix containing one or more columns of weights
#'
#' @return Returns a weighted mean for each combination of \code{y} by \code{x}.
#' @export
posterior_means = function(y, w) {
  y = as.matrix(y)
  w = as.matrix(w)
  as.tibble(t(y) %*% (w /colSums(w)))
}

### ----------------------------------------------------------



# Takes a set of integer weights representing a synthetic 
# population, expands them into the full population
#' Expand synthetic population weights
#'
#' @param wts Set of integer weights to expand
#'
#' @return Returns a vector where the index of each unit in the reference
#' sample is repeated \code{w} times where \code{w} is its weight.
#' @export
expand_sp_weights = function(wts) {
  wts %>% imap(function(wt, idx) {
    rep(idx, wt)
  }) %>%
    unlist()  
}

### ----------------------------------------------------------

#' Calculate odds of a proportion
#'
#' @param p A proportion
#'
#' @return Returns the odds of p as \code{p/(1-p)}
#' @export
odds = function(p) {
  p/(1-p)
}


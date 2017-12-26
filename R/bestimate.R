
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
#' @examples None yet
bestimate = function(samp,
                     ref,
                     y_var_names,
                     x_var_names,
                     sp_wts,
                     method = c("propensity", "prediction", "dr"),
                     posterior_draws = 100,
                     ... # Additional options passed to BART
                     ) {
  
  # Check to make sure that the X vars exist and are compatible across
  # samp and ref
  # Check to make sure y_vars exist 
  
  
  res = list()

  # Fit prediction model for each y_variable specified
  if ("prediction" %in% method) {
    res[["pred"]] = y_var_names %>%
      set_names(.) %>%
      map(~pbart2(x.train = select(samp, x_var_names), 
                  y.train = samp[[.]],
                  x.test = select(ref, x_var_names),
                  ...
                  )) %>%
      map("yhat.test") %>% # Extract the posterior for the ref sample
      map(pnorm) %>%       # Convert to probabilities
      map(~posterior_means(t(.), sp_wts)) %>%
      bind_rows(.id="y_var")
  }
  
  # If propensity or dr:
  # Get balanced subsample of the reference sample proportional to the weights
  
  
  
  # Fit propensity model to balanced samples
  # If prediction:
  # Fit a prediction model for each outcome variable
  # If DR
  # Append the mean propensity weight to the reference sample
  # Fit a DR model for each outcome variable/synthetic pop
  
  return(res)
}

# Utility function to multiply a matrix of 
# y values by a matrix of weights and return
# the weighted posterior distribution for y
posterior_means = function(y, w) {
  y = as.matrix(y)
  w = as.matrix(w)
  as.tibble(t(y) %*% (w /colSums(w)))
}



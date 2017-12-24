#' Title
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
#' @examples
bestimate = function(samp,
                     ref,
                     y_var_names,
                     x_var_names,
                     sp_wts=NULL,
                     method = c("propensity", "prediction", "dr"),
                     posterior_draws = 100) {
  
  # Check to make sure that the X vars exist and are compatible across
  # samp and ref
  # Check to make sure y_vars exist 
  
  # If propensity or dr:
  # Get balanced subsample of the reference sample proportional to the weights
  # Fit propensity model to balanced samples
  # If prediction:
  # Fit a prediction model for each outcome variable
  # If DR
  # Append the mean propensity weight to the reference sample
  # Fit a DR model for each outcome variable/synthetic pop
  
  
  
}
bestimate = function(samp,
                     ref,
                     y_var_names,
                     x_var_names,
                     replicate_weights=NULL,
                     method = c("propensity", "prediction", "dr"),
                     posterior_draws = 100) {
  
  # TODO Check to make sure that the X vars exist and are compatible across
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
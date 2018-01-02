#' @param y A data frame or matrix containing one or more columns of Y variables
#' @param w A data frame or matrix containing one or more columns of weights
#'
#' @return Returns a weighted mean for each combination of \code{y} by \code{x}.
#' @export

weighted_mean_posterior = function(y, w) {
  map2_dbl(as_tibble(y), as_tibble(w), weighted.mean)
}


ref_mean_components = function(y_pos, wts, name) {
  if (is.data.frame(y_pos)) {
    y_pos = list(y_pos)
  }
  wts %>%
    set_names(~sprintf("%s_%s", name, .)) %>%
    map(
      ~pmap(list(y = y_pos, w = .x), weighted_mean_posterior)
  ) %>% transpose() %>%
    map(bind_cols) 
    
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

#' Create subsample from a data.frame containing synthetic populations
#'
#' @param synth_pops \code{data.frame} containing synthetic populations
#' @param size Size of the subsamples
#' @param as_weights Indicates that the synthetic populations are 
#' represented as weights. If FALSE, synthetic populations are 
#' indices of records in the reference dataset.
#' @param seed Random seed
#'
#' @return A list containing subsamples for each synthetic population
#' @export
#'
#' @examples
sp_subsample = function(synth_pops, size, as_weights = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (as_weights) {
    synth_pops = map(synth_pops, expand_sp_weights)
  }
    map(synth_pops, ~sample(.x, size=size, replace=FALSE))
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

prop2weight = function(p, p_denom = TRUE) {
  if (p_denom) {
    w = odds(1-p)
  } else {
    w = odds(p)
  }
  w / sum(w)
}





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


#' Convert propensity to weight
#'
#' Takes a vector of propensity scores for a sample and
#' converts them into propensity weights using their
#' odds. By default, it assumes that the propensities
#' represent the probability than an observation comes from
#' the sample rather than the reference dataset. If the propensity
#' is the probability of being from the reference dataset, set
#' \code{p_denom = FALSE}.
#' 
#'
#' @param p A vector of propensity scores
#' @param p_denom If \code{TRUE}, indicates that the weights should
#' be calculated as \code{wt = (1 - p)/p}. Otherwise \code{wt = p/(1 - p)}
#'
#' @return A vector of propensity weights
#' @export
#'
#' @examples
prop2weight = function(p, p_denom = TRUE) {
  if (p_denom) {
    w = odds(1-p)
  } else {
    w = odds(p)
  }
  length(w) * w / sum(w)
}





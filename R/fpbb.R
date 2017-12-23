
#' @useDynLib bestimate
#' @importFrom Rcpp sourceCpp
NULL

#' Create synthetic populations with FPBB.
#'
#' @param weights Vector of weights for each unit in the reference sample. 
#'   These determine the initial probabilities of selection for the weighted
#'   Polya urn scheme used to generate each synthetic population.
#' @param L Number of synthetic populations to create.
#' @param N Size of the synthetic populations. In theory this is the size of the 
#'   full population but in practice only needs to be 10-20 times as large as
#'   the reference sample.
#' @param return_weights Logical flag indicating whether the synthetic
#'   populations should be returned as weights or indices. 
#'   \enumerate{
#'     \item If \code{TRUE}, output will contain the number of times each unit 
#'     in the reference sample appears in the synthetic population. 
#'     \item If \code{FALSE} output will be a vector containing the indices of 
#'     each unit in the order they were selected. This permits the use of 
#'     synthetic populations with procedure that cannot handle weights.
#'     }
#'
#' @return Returns a data frame with one column for each synthetic population.
#'   If \code{return_weights = TRUE} each column will be a vector containing
#'   the number of times a unit appears in the synthetic population. If 
#'   \code{return_weights = FALSE} each column will be a vector containing
#'   the row number of each unit in the reference sample in the order it was
#'   selected.
#' @export
#'
#' @examples
fpbb_synth_pops = function(weights, L = 2, 
                           N = length(weights) * 2, 
                           return_weights = TRUE) {
  
  n_ref = length(weights)
  
  # Rescale weights so that they sum to their length
  weights <- weights * n_ref/sum(weights)
  
  # Number of polya draws to create a synthetic population of size N
  num_draws = N - n_ref
  
  cat("Creating synthetic populations:\n")
  
  synth_pops = rerun(L, wtd_polya_sample_cpp(weights, num_draws)) %>%
    map(~c(., seq_len(n_ref)))
  
  return(synth_pops)
    
}  
  

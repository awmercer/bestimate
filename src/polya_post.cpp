
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector wtd_polya_sample_cpp(NumericVector wts, int size) {
    
    // Cumulative sum of weights
    NumericVector cum_wts = cumsum(wts);
    NumericVector::iterator wts_begin = cum_wts.begin();
    NumericVector::iterator wts_end = cum_wts.end();
    NumericVector::iterator wts_highest = wts_end - 1;
    
    // Random draws to be used to select sample
    NumericVector random_nums = runif(size);
    NumericVector::iterator random_nums_begin = random_nums.begin();
    NumericVector::iterator random_nums_end = random_nums.end();
    NumericVector::iterator rand;
    
    // Vector to hold the indices of the units 
    // that are drawn
    IntegerVector draws(size);
    IntegerVector::iterator draws_begin = draws.begin();
    IntegerVector::iterator draw;
    
    // This loop draws one unit according to the weighted
    // polya process 
    for(rand = random_nums_begin, draw = draws_begin; 
        rand != random_nums_end; 
        ++rand, ++draw) {
      
      // Make sure user hasn't interrupted
      Rcpp::checkUserInterrupt();
      

      
      // Get random position in cumulative weights
      int draw_value = *rand * *wts_highest;
      NumericVector::iterator pos;
      pos = std::upper_bound(wts_begin, wts_end, draw_value);

      // Set draw equal to the location of the draw in the cumulative weights
      // Add 1 to the value so that it is consistent with R indexing
      // when the vector is passed back to R.
      *draw = std::distance(wts_begin, pos) + 1;
      
      // Increment the cumulative weights by 1 from pos to end
      for(NumericVector::iterator i = pos; i != wts_end; ++i) {
        ++*i;
      }
    }

    return draws;
}

/*** R
wtd_polya_sample_cpp(size=5, wts=c(4,6,10))
  */



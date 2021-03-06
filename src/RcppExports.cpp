// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// wtd_polya_sample_cpp
IntegerVector wtd_polya_sample_cpp(NumericVector wts, int size);
RcppExport SEXP _bestimate_wtd_polya_sample_cpp(SEXP wtsSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(wtd_polya_sample_cpp(wts, size));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bestimate_wtd_polya_sample_cpp", (DL_FUNC) &_bestimate_wtd_polya_sample_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_bestimate(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

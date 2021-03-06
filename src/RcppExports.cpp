// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// SEB
List SEB(arma::mat X, arma::uvec nints, const bool intervals, const std::string order);
RcppExport SEXP _SEB_SEB(SEXP XSEXP, SEXP nintsSEXP, SEXP intervalsSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type nints(nintsSEXP);
    Rcpp::traits::input_parameter< const bool >::type intervals(intervalsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(SEB(X, nints, intervals, order));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SEB_SEB", (DL_FUNC) &_SEB_SEB, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SEB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

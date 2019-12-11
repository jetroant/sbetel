// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// test_rcpp
double test_rcpp(double num);
RcppExport SEXP _sbetel_test_rcpp(SEXP numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type num(numSEXP);
    rcpp_result_gen = Rcpp::wrap(test_rcpp(num));
    return rcpp_result_gen;
END_RCPP
}
// etel_rcpp
double etel_rcpp(NumericVector th, Function g, int p, NumericMatrix y, int itermax, int bw, double lambda, int td);
RcppExport SEXP _sbetel_etel_rcpp(SEXP thSEXP, SEXP gSEXP, SEXP pSEXP, SEXP ySEXP, SEXP itermaxSEXP, SEXP bwSEXP, SEXP lambdaSEXP, SEXP tdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type th(thSEXP);
    Rcpp::traits::input_parameter< Function >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< int >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type td(tdSEXP);
    rcpp_result_gen = Rcpp::wrap(etel_rcpp(th, g, p, y, itermax, bw, lambda, td));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sbetel_test_rcpp", (DL_FUNC) &_sbetel_test_rcpp, 1},
    {"_sbetel_etel_rcpp", (DL_FUNC) &_sbetel_etel_rcpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_sbetel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
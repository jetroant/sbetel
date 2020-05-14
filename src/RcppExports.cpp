// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// etel_rcpp
double etel_rcpp(NumericVector th, Function g, NumericMatrix y, int bw, int td, int itermax, List args);
RcppExport SEXP _sbetel_etel_rcpp(SEXP thSEXP, SEXP gSEXP, SEXP ySEXP, SEXP bwSEXP, SEXP tdSEXP, SEXP itermaxSEXP, SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type th(thSEXP);
    Rcpp::traits::input_parameter< Function >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type td(tdSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(etel_rcpp(th, g, y, bw, td, itermax, args));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sbetel_etel_rcpp", (DL_FUNC) &_sbetel_etel_rcpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_sbetel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

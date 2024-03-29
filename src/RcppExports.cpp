// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bhnoneq_LL
double bhnoneq_LL(NumericVector stpar, NumericVector year, NumericVector Lbar, NumericVector ss, double Linf, double K, double Lc, int nbreaks);
RcppExport SEXP _DLMtool_bhnoneq_LL(SEXP stparSEXP, SEXP yearSEXP, SEXP LbarSEXP, SEXP ssSEXP, SEXP LinfSEXP, SEXP KSEXP, SEXP LcSEXP, SEXP nbreaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type stpar(stparSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type year(yearSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lbar(LbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< double >::type Linf(LinfSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type Lc(LcSEXP);
    Rcpp::traits::input_parameter< int >::type nbreaks(nbreaksSEXP);
    rcpp_result_gen = Rcpp::wrap(bhnoneq_LL(stpar, year, Lbar, ss, Linf, K, Lc, nbreaks));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DLMtool_bhnoneq_LL", (DL_FUNC) &_DLMtool_bhnoneq_LL, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_DLMtool(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

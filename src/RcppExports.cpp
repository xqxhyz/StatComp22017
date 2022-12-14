// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gibbsC
NumericMatrix gibbsC(double mu1, double mu2, double sigma1, double sigma2, double rho, int N, int burn);
RcppExport SEXP _StatComp22017_gibbsC(SEXP mu1SEXP, SEXP mu2SEXP, SEXP sigma1SEXP, SEXP sigma2SEXP, SEXP rhoSEXP, SEXP NSEXP, SEXP burnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< double >::type mu2(mu2SEXP);
    Rcpp::traits::input_parameter< double >::type sigma1(sigma1SEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsC(mu1, mu2, sigma1, sigma2, rho, N, burn));
    return rcpp_result_gen;
END_RCPP
}
// til
double til(double LYM, double STR, double TUM, double TIL3, double TIL7, double TIL8);
RcppExport SEXP _StatComp22017_til(SEXP LYMSEXP, SEXP STRSEXP, SEXP TUMSEXP, SEXP TIL3SEXP, SEXP TIL7SEXP, SEXP TIL8SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type LYM(LYMSEXP);
    Rcpp::traits::input_parameter< double >::type STR(STRSEXP);
    Rcpp::traits::input_parameter< double >::type TUM(TUMSEXP);
    Rcpp::traits::input_parameter< double >::type TIL3(TIL3SEXP);
    Rcpp::traits::input_parameter< double >::type TIL7(TIL7SEXP);
    Rcpp::traits::input_parameter< double >::type TIL8(TIL8SEXP);
    rcpp_result_gen = Rcpp::wrap(til(LYM, STR, TUM, TIL3, TIL7, TIL8));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp22017_gibbsC", (DL_FUNC) &_StatComp22017_gibbsC, 7},
    {"_StatComp22017_til", (DL_FUNC) &_StatComp22017_til, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp22017(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

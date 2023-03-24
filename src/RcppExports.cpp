// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// TS_cont
Rcpp::NumericVector TS_cont(Rcpp::NumericVector x, Rcpp::NumericVector Fx, Rcpp::NumericVector param, Rcpp::Function qnull, Rcpp::CharacterVector doMethod);
RcppExport SEXP _Rgof_TS_cont(SEXP xSEXP, SEXP FxSEXP, SEXP paramSEXP, SEXP qnullSEXP, SEXP doMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Fx(FxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type qnull(qnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type doMethod(doMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(TS_cont(x, Fx, param, qnull, doMethod));
    return rcpp_result_gen;
END_RCPP
}
// TS_disc
NumericVector TS_disc(IntegerVector x, NumericVector p, NumericMatrix nm, NumericVector vals, CharacterVector doMethod);
RcppExport SEXP _Rgof_TS_disc(SEXP xSEXP, SEXP pSEXP, SEXP nmSEXP, SEXP valsSEXP, SEXP doMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nm(nmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type doMethod(doMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(TS_disc(x, p, nm, vals, doMethod));
    return rcpp_result_gen;
END_RCPP
}
// bincounter_cpp
Rcpp::IntegerVector bincounter_cpp(Rcpp::NumericVector x, Rcpp::NumericVector bins);
RcppExport SEXP _Rgof_bincounter_cpp(SEXP xSEXP, SEXP binsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type bins(binsSEXP);
    rcpp_result_gen = Rcpp::wrap(bincounter_cpp(x, bins));
    return rcpp_result_gen;
END_RCPP
}
// binner_cont
Rcpp::NumericVector binner_cont(Rcpp::NumericVector x, Rcpp::Function pnull, Rcpp::NumericVector param, int k, int which, Rcpp::NumericVector Range);
RcppExport SEXP _Rgof_binner_cont(SEXP xSEXP, SEXP pnullSEXP, SEXP paramSEXP, SEXP kSEXP, SEXP whichSEXP, SEXP RangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type which(whichSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Range(RangeSEXP);
    rcpp_result_gen = Rcpp::wrap(binner_cont(x, pnull, param, k, which, Range));
    return rcpp_result_gen;
END_RCPP
}
// binner_disc
Rcpp::IntegerVector binner_disc(Rcpp::IntegerVector x, Rcpp::NumericVector p, int k);
RcppExport SEXP _Rgof_binner_disc(SEXP xSEXP, SEXP pSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(binner_disc(x, p, k));
    return rcpp_result_gen;
END_RCPP
}
// chi_stat_cont
double chi_stat_cont(Rcpp::NumericVector param, Rcpp::NumericVector x, Rcpp::Function pnull, Rcpp::NumericVector bins, std::string formula, double rate);
RcppExport SEXP _Rgof_chi_stat_cont(SEXP paramSEXP, SEXP xSEXP, SEXP pnullSEXP, SEXP binsSEXP, SEXP formulaSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< std::string >::type formula(formulaSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(chi_stat_cont(param, x, pnull, bins, formula, rate));
    return rcpp_result_gen;
END_RCPP
}
// chi_stat_disc
double chi_stat_disc(Rcpp::NumericVector param, Rcpp::IntegerVector x, Rcpp::Function pnull, Rcpp::IntegerVector bins, std::string formula, double rate);
RcppExport SEXP _Rgof_chi_stat_disc(SEXP paramSEXP, SEXP xSEXP, SEXP pnullSEXP, SEXP binsSEXP, SEXP formulaSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< std::string >::type formula(formulaSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(chi_stat_disc(param, x, pnull, bins, formula, rate));
    return rcpp_result_gen;
END_RCPP
}
// chi_test_cont
Rcpp::NumericMatrix chi_test_cont(Rcpp::NumericVector x, Rcpp::Function pnull, Rcpp::NumericVector param, std::string formula, double rate, Rcpp::IntegerVector nbins, Rcpp::NumericVector Range, int Minimize);
RcppExport SEXP _Rgof_chi_test_cont(SEXP xSEXP, SEXP pnullSEXP, SEXP paramSEXP, SEXP formulaSEXP, SEXP rateSEXP, SEXP nbinsSEXP, SEXP RangeSEXP, SEXP MinimizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< std::string >::type formula(formulaSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nbins(nbinsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Range(RangeSEXP);
    Rcpp::traits::input_parameter< int >::type Minimize(MinimizeSEXP);
    rcpp_result_gen = Rcpp::wrap(chi_test_cont(x, pnull, param, formula, rate, nbins, Range, Minimize));
    return rcpp_result_gen;
END_RCPP
}
// chi_test_disc
Rcpp::NumericMatrix chi_test_disc(Rcpp::IntegerVector x, Rcpp::Function pnull, Rcpp::NumericVector param, Rcpp::IntegerVector nbins, std::string formula, double rate, int Minimize);
RcppExport SEXP _Rgof_chi_test_disc(SEXP xSEXP, SEXP pnullSEXP, SEXP paramSEXP, SEXP nbinsSEXP, SEXP formulaSEXP, SEXP rateSEXP, SEXP MinimizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nbins(nbinsSEXP);
    Rcpp::traits::input_parameter< std::string >::type formula(formulaSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< int >::type Minimize(MinimizeSEXP);
    rcpp_result_gen = Rcpp::wrap(chi_test_disc(x, pnull, param, nbins, formula, rate, Minimize));
    return rcpp_result_gen;
END_RCPP
}
// gof_cont
Rcpp::NumericMatrix gof_cont(Rcpp::NumericVector x, Rcpp::Function pnull, Rcpp::Function phat, Rcpp::Function rnull, Rcpp::Function qnull, int B, Rcpp::CharacterVector doMethod);
RcppExport SEXP _Rgof_gof_cont(SEXP xSEXP, SEXP pnullSEXP, SEXP phatSEXP, SEXP rnullSEXP, SEXP qnullSEXP, SEXP BSEXP, SEXP doMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type rnull(rnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type qnull(qnullSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type doMethod(doMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(gof_cont(x, pnull, phat, rnull, qnull, B, doMethod));
    return rcpp_result_gen;
END_RCPP
}
// gof_disc
NumericMatrix gof_disc(Rcpp::IntegerVector x, Rcpp::Function pnull, Rcpp::NumericVector vals, Rcpp::Function rnull, Rcpp::Function phat, double rate, int B, Rcpp::CharacterVector doMethod);
RcppExport SEXP _Rgof_gof_disc(SEXP xSEXP, SEXP pnullSEXP, SEXP valsSEXP, SEXP rnullSEXP, SEXP phatSEXP, SEXP rateSEXP, SEXP BSEXP, SEXP doMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type rnull(rnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type doMethod(doMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(gof_disc(x, pnull, vals, rnull, phat, rate, B, doMethod));
    return rcpp_result_gen;
END_RCPP
}
// nm_calc
NumericMatrix nm_calc(int n);
RcppExport SEXP _Rgof_nm_calc(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(nm_calc(n));
    return rcpp_result_gen;
END_RCPP
}
// power_cont
Rcpp::NumericMatrix power_cont(Rcpp::Function pnull, Rcpp::Function phat, Rcpp::Function rnull, Rcpp::Function ralt, Rcpp::NumericVector param_alt, Rcpp::Function qnull, Rcpp::IntegerVector nbins, double rate, Rcpp::NumericVector Range, Rcpp::IntegerVector B, const double alpha);
RcppExport SEXP _Rgof_power_cont(SEXP pnullSEXP, SEXP phatSEXP, SEXP rnullSEXP, SEXP raltSEXP, SEXP param_altSEXP, SEXP qnullSEXP, SEXP nbinsSEXP, SEXP rateSEXP, SEXP RangeSEXP, SEXP BSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type rnull(rnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type ralt(raltSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param_alt(param_altSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type qnull(qnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nbins(nbinsSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Range(RangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(power_cont(pnull, phat, rnull, ralt, param_alt, qnull, nbins, rate, Range, B, alpha));
    return rcpp_result_gen;
END_RCPP
}
// power_disc
Rcpp::NumericMatrix power_disc(Rcpp::Function pnull, Rcpp::Function phat, Rcpp::Function rnull, Rcpp::Function ralt, Rcpp::NumericVector param_alt, Rcpp::NumericVector vals, Rcpp::IntegerVector nbins, double rate, Rcpp::IntegerVector B, const double alpha);
RcppExport SEXP _Rgof_power_disc(SEXP pnullSEXP, SEXP phatSEXP, SEXP rnullSEXP, SEXP raltSEXP, SEXP param_altSEXP, SEXP valsSEXP, SEXP nbinsSEXP, SEXP rateSEXP, SEXP BSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type rnull(rnullSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type ralt(raltSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param_alt(param_altSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nbins(nbinsSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(power_disc(pnull, phat, rnull, ralt, param_alt, vals, nbins, rate, B, alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rgof_TS_cont", (DL_FUNC) &_Rgof_TS_cont, 5},
    {"_Rgof_TS_disc", (DL_FUNC) &_Rgof_TS_disc, 5},
    {"_Rgof_bincounter_cpp", (DL_FUNC) &_Rgof_bincounter_cpp, 2},
    {"_Rgof_binner_cont", (DL_FUNC) &_Rgof_binner_cont, 6},
    {"_Rgof_binner_disc", (DL_FUNC) &_Rgof_binner_disc, 3},
    {"_Rgof_chi_stat_cont", (DL_FUNC) &_Rgof_chi_stat_cont, 6},
    {"_Rgof_chi_stat_disc", (DL_FUNC) &_Rgof_chi_stat_disc, 6},
    {"_Rgof_chi_test_cont", (DL_FUNC) &_Rgof_chi_test_cont, 8},
    {"_Rgof_chi_test_disc", (DL_FUNC) &_Rgof_chi_test_disc, 7},
    {"_Rgof_gof_cont", (DL_FUNC) &_Rgof_gof_cont, 7},
    {"_Rgof_gof_disc", (DL_FUNC) &_Rgof_gof_disc, 8},
    {"_Rgof_nm_calc", (DL_FUNC) &_Rgof_nm_calc, 1},
    {"_Rgof_power_cont", (DL_FUNC) &_Rgof_power_cont, 11},
    {"_Rgof_power_disc", (DL_FUNC) &_Rgof_power_disc, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rgof(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

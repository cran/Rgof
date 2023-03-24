#ifndef CHI_STAT_CONT_H
#define CHI_STAT_CONT_H

#include <Rcpp.h>
double chi_stat_cont(
      Rcpp::NumericVector param,
      Rcpp::NumericVector x,
      Rcpp::Function pnull,
      Rcpp::NumericVector bins,
      std::string formula,
      double rate);

#endif

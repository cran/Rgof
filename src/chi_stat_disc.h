#ifndef CHI_STAT_DISC_H
#define CHI_STAT_DISC_H

#include <Rcpp.h>
double chi_stat_disc(
      Rcpp::NumericVector pram,
      Rcpp::IntegerVector x,
      Rcpp::Function pnull,
      Rcpp::IntegerVector bins,
      std::string formula,
      double rate);

#endif

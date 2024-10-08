#ifndef WBINCOUNTER_H
#define WBINCOUNTER_H

#include <Rcpp.h>
Rcpp::NumericVector wbincounter(
      Rcpp::NumericVector x, 
      Rcpp::NumericVector bins,
      Rcpp::NumericVector w);

#endif

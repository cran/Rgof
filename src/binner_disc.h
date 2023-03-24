#ifndef BINNER_DISC_H
#define BINNER_DISC_H

#include <Rcpp.h>
Rcpp::IntegerVector binner_disc(
      Rcpp::IntegerVector x, 
      Rcpp::NumericVector p, 
      int k=10);

#endif

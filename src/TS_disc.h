#ifndef TS_DISC_H
#define TS_DISC_H

#include <Rcpp.h>
Rcpp::NumericVector TS_disc(Rcpp::IntegerVector x, 
                            Rcpp::NumericVector p,
                            Rcpp::NumericMatrix nm,
                            Rcpp::NumericVector vals,
                            Rcpp::CharacterVector doMethod
                      );

#endif

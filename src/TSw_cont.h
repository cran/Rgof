#ifndef TSw_CONT_H
#define TSw_CONT_H

#include <Rcpp.h>
Rcpp::NumericVector TSw_cont(Rcpp::NumericVector x, 
                            Rcpp::Function pnull,
                            Rcpp::NumericVector w
                      );

#endif

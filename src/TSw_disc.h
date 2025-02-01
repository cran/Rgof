#ifndef TSW_DISC_H
#define TSW_DISC_H

#include <Rcpp.h>
Rcpp::NumericVector TSw_disc(Rcpp::IntegerVector x, 
                            Rcpp::Function pnull,
                            Rcpp::NumericVector param,
                            Rcpp::NumericVector vals,
                            Rcpp::NumericVector w);

#endif

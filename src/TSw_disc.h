#ifndef TSW_DISC_H
#define TSW_DISC_H

#include <Rcpp.h>
Rcpp::NumericVector TSw_disc(Rcpp::IntegerVector x, 
                            Rcpp::NumericVector Fx,
                            Rcpp::NumericVector vals ,
                            Rcpp::NumericVector w);

#endif

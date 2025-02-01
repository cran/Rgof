#ifndef TS_DISC_H
#define TS_DISC_H

#include <Rcpp.h>
Rcpp::NumericVector TS_disc(Rcpp::IntegerVector x, 
                            Rcpp::Function pnull,
                            Rcpp::NumericVector param,
                            Rcpp::NumericVector vals
                      );

#endif

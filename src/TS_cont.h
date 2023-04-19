#ifndef TS_CONT_H
#define TS_CONT_H

#include <Rcpp.h>
Rcpp::NumericVector TS_cont(Rcpp::NumericVector x, 
                            Rcpp::NumericVector Fx,
                            Rcpp::NumericVector param, 
                            Rcpp::Function qnull
                      );

#endif

#ifndef GOF_DISC_H
#define GOF_DISC_H

#include <Rcpp.h>
Rcpp::NumericMatrix gof_disc(Rcpp::IntegerVector x, 
                       Rcpp::Function pnull, 
                       Rcpp::Function rnull, 
                       Rcpp::NumericVector vals,                        
                       Rcpp::Function phat,
                       Rcpp::Function TS,
                       int NewTS,
                       Rcpp::List TSextra,
                       double rate=0, 
                       int B=5000);

#endif

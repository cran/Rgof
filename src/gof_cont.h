#ifndef GOF_CONT_H
#define GOF_CONT_H

#include <Rcpp.h>
Rcpp::NumericMatrix gof_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::Function rnull, 
        Rcpp::Function qnull,
        Rcpp::Function phat, 
        Rcpp::Function TS_cont_alt,        
        int B=5000
        );

#endif

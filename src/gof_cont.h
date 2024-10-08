#ifndef GOF_CONT_H
#define GOF_CONT_H

#include <Rcpp.h>
Rcpp::NumericMatrix gof_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::Function rnull, 
        Rcpp::Function qnull,
        Rcpp::Function w,
        Rcpp::Function phat, 
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,        
        int B=5000
        );

#endif

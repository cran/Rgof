#ifndef GOF_CONT_H
#define GOF_CONT_H

#include <Rcpp.h>
Rcpp::NumericMatrix gof_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::Function phat,
        Rcpp::Function rnull, 
        Rcpp::Function qnull, 
        int B=5000,
        Rcpp::CharacterVector doMethod=Rcpp::CharacterVector::create("KS", "Kuiper", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1")
        );

#endif

#ifndef GOF_DISC_H
#define GOF_DISC_H

#include <Rcpp.h>
Rcpp::NumericMatrix gof_disc(Rcpp::IntegerVector x, 
                       Rcpp::Function pnull, 
                       Rcpp::NumericVector vals,                        
                       Rcpp::Function rnull, 
                       Rcpp::Function phat,
                       double rate=0, 
                       int B=5000,
             Rcpp::CharacterVector doMethod=Rcpp::CharacterVector::create("KS", "Kuiper", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1"));

#endif

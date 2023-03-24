#ifndef CHI_TEST_CONT_H
#define CHI_TEST_CONT_H

#include <Rcpp.h>
Rcpp::NumericMatrix chi_test_cont (
    Rcpp::NumericVector x, 
    Rcpp::Function pnull, 
    Rcpp::NumericVector param,
    std::string formula,
    double rate=0.0,
    Rcpp::IntegerVector nbins= Rcpp::IntegerVector::create(100, 10),
    Rcpp::NumericVector Range=Rcpp::NumericVector::create(-99, -99),
    int Minimize=0);

#endif

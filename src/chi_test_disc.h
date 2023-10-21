#ifndef CHI_TEST_DISC_H
#define CHI_TEST_DISC_H

#include <Rcpp.h>
Rcpp::NumericMatrix chi_test_disc (
    Rcpp::IntegerVector x, 
    Rcpp::Function pnull, 
    Rcpp::NumericVector param, 
    Rcpp::IntegerVector nbins= Rcpp::IntegerVector::create(100, 10),
    std::string formula="Pearson",
    double rate=0.0,
    int Minimize=0,
    double minexpcount=2.0
  );

#endif

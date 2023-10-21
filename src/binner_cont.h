#ifndef BINNER_CONT_H
#define BINNER_CONT_H

#include <Rcpp.h>
Rcpp::NumericVector binner_cont(
      Rcpp::NumericVector x, 
      Rcpp::Function pnull,
      Rcpp::NumericVector param,
      int k=10,
      int which=1,
      Rcpp::NumericVector Range=Rcpp::NumericVector::create(-99, -99),
      double minexpcount=2.0);

#endif

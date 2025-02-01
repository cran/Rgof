#include <Rcpp.h>
#include "gof_disc.h"
#include "TS_disc.h"

using namespace Rcpp;

//' find power of gof tests for discrete data
//' 
//' @param pnull R function (cdf)
//' @param rnull R function (generate data under null hypothesis)
//' @param vals vector of values of discrete random variable
//' @param ralt  R function to generate data under alternative
//' @param param_alt parameters of function ralt
//' @param phat  function to estimate parameters from the data
//' @param TS function to calculate test statistics
//' @param typeTS type of test statistic
//' @param TSextra list passed to TS, if desired
//' @param rate =0, rate of sample size, if random
//' @param B  =c(1000, 1000) Number of simulation runs for power and null distribution
//' @param alpha =0.05, type I error of test
//' @keywords internal
//' @return A matrix of powers
// [[Rcpp::export]]
Rcpp::NumericMatrix power_disc(
        Rcpp::Function pnull, 
        Rcpp::Function rnull, 
        Rcpp::NumericVector vals,         
        Rcpp::Function ralt, 
        Rcpp::NumericVector param_alt,
        Rcpp::Function phat, 
        Rcpp::Function TS,
        int typeTS, 
        Rcpp::List TSextra,
        double rate=0.0,
        Rcpp::IntegerVector B=Rcpp::IntegerVector::create(1000, 1000), 
        const double alpha=0.05) {
  
  int  i, j, k, np=param_alt.size();
  IntegerVector x=ralt(param_alt(0));
  NumericVector TS_data;
  
  if(typeTS<=1) TS_data=TS(x, pnull, phat(x), vals); 
  if(typeTS==2) TS_data=TS(x, pnull, phat(x), vals, TSextra); 
  int const nummethods=TS_data.size();
  NumericMatrix out(np, nummethods);
  for(i=0;i<B(0);++i) {
    for(j=0;j<np;++j) {
      x=ralt(param_alt[j]); 
      NumericMatrix tmp = gof_disc(x, pnull, rnull, vals, phat, 
              TS, typeTS, TSextra, rate, B[1]);
      for(k=0;k<nummethods;++k) 
         if(tmp(1,k)<alpha)  out(j, k) = out(j, k)+1;
      }
  } 
  return out/B(0);
}

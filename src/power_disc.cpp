#include <Rcpp.h>
#include "gof_disc.h"
#include "nm_calc.h"
#include "TS_disc.h"
#include "chi_test_disc.h"

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
//' @param nbins =c(100,10) number of bins to use 
//' @param rate rate of Poisson if sample size is random, 0 otherwise
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
        Rcpp::IntegerVector nbins=Rcpp::IntegerVector::create(100, 10), 
        double rate=0.0,
        Rcpp::IntegerVector B=Rcpp::IntegerVector::create(1000, 1000), 
        const double alpha=0.05) {
  
  int  i, j, k, np=param_alt.size();
  IntegerVector x=ralt(param_alt(0));
  NumericMatrix nm=nm_calc(sum(x));
  NumericVector Fx(x.size());
  for(i=0;i<x.size();++i) Fx(i)=(1.0+i)/(1.0+i);
  NumericVector TS_data=TS(x, Fx, nm, vals);  
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector tsMethods=TS_data.names();
  Rcpp::CharacterVector chiMethods= CharacterVector::create( 
    "chi large Pearson", "chi small Pearson",  "chi large LR", "chi small LR");  
  Rcpp::CharacterVector doMethod(nummethods+4);
  for(i=0;i<nummethods;++i) doMethod[i]=tsMethods[i];
  for(i=0;i<4;++i) doMethod[i+nummethods]=chiMethods[i];  

  NumericMatrix out(np, nummethods+4);
  for(i=0;i<B(0);++i) {
    for(j=0;j<np;++j) {
      IntegerVector x=ralt(param_alt[j]); 
      NumericMatrix tmp1 = gof_disc(x, pnull, rnull, vals, phat, TS, rate, B(1));
      for(k=0;k<8;++k) 
        if(tmp1(1,k)<alpha)  out(j, k) = out(j, k)+1;
      NumericVector param=phat(x);           
      NumericMatrix tmp2 = chi_test_disc(x, pnull, param,  
                nbins, "Pearson", rate, 1);
      for(k=0;k<2;++k) 
        if(tmp2(k, 2)<alpha) out(j, nummethods+k) = out(j, nummethods+k)+1; 
      tmp2 = chi_test_disc(x, pnull, param, nbins, "LR", rate, 1);
      for(k=0;k<2;++k) 
          if(tmp2(k, 2)<alpha) out(j, nummethods+2+k) = out(j, nummethods+2+k)+1;   
  }    
} 
  return out/B(0);
}

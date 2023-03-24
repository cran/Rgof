#include <Rcpp.h>
#include "gof_disc.h"
#include "chi_test_disc.h"

using namespace Rcpp;

//' find power of gof tests for discrete data
//' 
//' @param pnull R function (cdf)
//' @param phat  function to estimate parameters from the data
//' @param rnull R function (generate data under null hypothesis)
//' @param ralt  R function to generate data under alternative
//' @param param_alt parameters of function ralt
//' @param vals vector of values of discrete random variable
//' @param nbins =c(100,10) number of bins to use 
//' @param rate rate of Poisson if sample size is random, 0 otherwise
//' @param B  =c(1000, 1000) Number of simulation runs for power and null distribution
//' @param alpha =0.05, type I error of test 
//' @keywords internal
//' @return A matrix of powers
// [[Rcpp::export]]
Rcpp::NumericMatrix power_disc(
        Rcpp::Function pnull, 
        Rcpp::Function phat,
        Rcpp::Function rnull, 
        Rcpp::Function ralt, 
        Rcpp::NumericVector param_alt,
        Rcpp::NumericVector vals, 
        Rcpp::IntegerVector nbins=Rcpp::IntegerVector::create(100, 10), 
        double rate=0.0,
        Rcpp::IntegerVector B=Rcpp::IntegerVector::create(1000, 1000), 
        const double alpha=0.05) {
  
  int  i, j, k, np=param_alt.size();
  CharacterVector doMethod= CharacterVector::create("KS", 
          "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", 
          "Wassp1", "chi large Pearson", "chi small Pearson",
          "chi large LR", "chi small LR");

  NumericMatrix out(np, 13);
  for(i=0;i<B(0);++i) {
    for(j=0;j<np;++j) {
      IntegerVector x=ralt(param_alt[j]); 
      NumericMatrix tmp1 = gof_disc(x, pnull, vals, rnull, phat, rate, B(1));
      for(k=0;k<9;++k) 
        if(tmp1(1,k)<alpha)  out(j, k) = out(j, k)+1;
      NumericVector param=phat(x);           
      NumericMatrix tmp2 = chi_test_disc(x, pnull, param,  
                nbins, "Pearson", rate, 1);
      for(k=0;k<2;++k) 
        if(tmp2(k, 2)<alpha) out(j, k+9) = out(j, k+9)+1; 
      tmp2 = chi_test_disc(x, pnull, param, nbins, "LR", rate, 1);
      for(k=0;k<2;++k) 
          if(tmp2(k, 2)<alpha) out(j, k+11) = out(j, k+11)+1;   
  }    
} 
  return out/B(0);
}

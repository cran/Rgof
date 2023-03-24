#include <Rcpp.h>
#include <string>
#include "gof_cont.h"
#include "chi_test_cont.h"

using namespace Rcpp;

//' find power of gof tests for continuous data
//' 
//' @param pnull R function (cdf)
//' @param phat  function to estimate parameters from the data
//' @param rnull R function (generate data under null hypothesis)
//' @param ralt  R function to generate data under alternative
//' @param param_alt parameters of ralt
//' @param qnull R function (quantiles under null hypothesis)
//' @param nbins =c(100,10) number of bins to use
//' @param rate rate of Poisson if sample size is random
//' @param Range =(-99999, 99999) limits of possible observations, if any
//' @param B  =c(1000, 1000) Number of simulation runs for power and null distribution
//' @param alpha =0.05, type I error of test 
//' @keywords internal
//' @return A matrix of powers
// [[Rcpp::export]]
Rcpp::NumericMatrix power_cont(
        Rcpp::Function pnull, 
        Rcpp::Function phat,
        Rcpp::Function rnull, 
        Rcpp::Function ralt, 
        Rcpp::NumericVector param_alt,
        Rcpp::Function qnull, 
        Rcpp::IntegerVector nbins=Rcpp::IntegerVector::create(100, 10), 
        double rate=0.0,
        Rcpp::NumericVector Range=Rcpp::NumericVector::create(-99999, 99999),
        Rcpp::IntegerVector B=Rcpp::IntegerVector::create(1000, 1000), 
        const double alpha=0.05) {
        
  int const nummethods=17;
  int  i, j, k, np=param_alt.size();
  CharacterVector doMethod= CharacterVector::create("KS", 
         "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1", 
         "EP large Pearson", "ES large Pearson", "EP small Pearson", "ES small Pearson",
         "EP large LR", "ES large LR", "EP small LR", "ES small LR");
  NumericMatrix out(np, nummethods);

  for(i=0;i<B(0);++i) {
     for(j=0;j<np;++j) {
         NumericVector x=ralt(param_alt[j]); 
         NumericMatrix tmp1 = gof_cont(x, pnull, phat, 
                   rnull,  qnull, B(1), doMethod);
         for(k=0;k<9;++k) 
            if(tmp1(1,k)<alpha) out(j, k) = out(j, k)+1;            
         NumericVector param=phat(x);     
         NumericMatrix tmp2 = chi_test_cont(x, pnull, 
                  param, "Pearson", rate, nbins, Range, 1);        
         for(k=0;k<4;++k) 
           if(tmp2(k, 2)<alpha) out(j, k+9) = out(j, k+9)+1;
         tmp2 = chi_test_cont(x, pnull, 
                  param, "LR", rate, nbins, Range, 1);        
         for(k=0;k<4;++k) 
             if(tmp2(k, 2)<alpha) out(j, k+13) = out(j, k+13)+1;           
         }
  }
  
  return out/B(0);
}

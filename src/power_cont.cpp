#include <Rcpp.h>
#include <string>
#include "gof_cont.h"
#include "TS_cont.h"
#include "chi_test_cont.h"

using namespace Rcpp;

//' find power of gof tests for continuous data
//' 
//' @param pnull R function (cdf)
//' @param rnull R function (generate data under null hypothesis)
//' @param qnull R function (quantiles under null hypothesis)
//' @param ralt  R function to generate data under alternative
//' @param param_alt parameters of ralt
//' @param phat  function to estimate parameters from the data
//' @param TS function to calculate test statistics
//' @param nbins =c(100,10) number of bins to use
//' @param rate rate of Poisson if sample size is random
//' @param Range =(-99999, 99999) limits of possible observations, if any
//' @param B  =c(1000, 1000) Number of simulation runs for power and null distribution
//' @param alpha =0.05, type I error of test 
//' @param minexpcount =2 minimal expected bin count required
//' @keywords internal
//' @return A matrix of powers
// [[Rcpp::export]]
Rcpp::NumericMatrix power_cont(
        Rcpp::Function pnull, 
        Rcpp::Function rnull,       
        Rcpp::Function qnull, 
        Rcpp::Function ralt, 
        Rcpp::NumericVector param_alt,
        Rcpp::Function phat,  
        Rcpp::Function TS,
        Rcpp::IntegerVector nbins=Rcpp::IntegerVector::create(100, 10), 
        double rate=0.0,
        Rcpp::NumericVector Range=Rcpp::NumericVector::create(-99999, 99999),
        Rcpp::IntegerVector B=Rcpp::IntegerVector::create(1000, 1000), 
        const double alpha=0.05,
        double minexpcount=2.0) {
  
  int  i, j, k, np=param_alt.size();
  NumericVector x=ralt(param_alt[0]);
  NumericVector Fx(x.size());
  for(i=0;i<x.length();++i) Fx(i)=(1.0+i)/(2.0+i);
  NumericVector TS_data=TS(x, Fx, 1.0, qnull);  
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector tsMethods=TS_data.names();
  Rcpp::CharacterVector chiMethods= CharacterVector::create("EP large Pearson", 
              "ES large Pearson", "EP small Pearson", "ES small Pearson",
              "EP large LR", "ES large LR", "EP small LR", "ES small LR");  
  Rcpp::CharacterVector doMethod(nummethods+8);
  for(i=0;i<nummethods;++i) doMethod[i]=tsMethods[i];
  for(i=0;i<8;++i) doMethod[i+nummethods]=chiMethods[i];
  NumericMatrix out(np, nummethods+8);
  colnames(out) = doMethod;
  for(i=0;i<B(0);++i) {
     for(j=0;j<np;++j) {
         NumericVector x=ralt(param_alt[j]); 
         NumericMatrix tmp1 = gof_cont(x, pnull,  
                   rnull,  qnull, phat, TS, B(1));
         for(k=0;k<nummethods;++k) 
            if(tmp1(1,k)<alpha) out(j, k) = out(j, k)+1;            
         NumericVector param=phat(x);     
         NumericMatrix tmp2 = chi_test_cont(x, pnull, 
                  param, "Pearson", rate, nbins, Range, 1, minexpcount);        
         for(k=0;k<4;++k) 
           if(tmp2(k, 2)<alpha) out(j, nummethods+k) = out(j, nummethods+k)+1;
         tmp2 = chi_test_cont(x, pnull, 
                  param, "LR", rate, nbins, Range, 1, minexpcount);        
         for(k=0;k<4;++k) 
             if(tmp2(k, 2)<alpha) out(j, nummethods+4+k) = out(j, nummethods+4+k)+1;           
         }
  }
  
  return out/B(0);
}

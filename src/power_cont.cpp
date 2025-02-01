#include <Rcpp.h>
#include <string>
#include "gof_cont.h"

using namespace Rcpp;

//' find power of gof tests for continuous data
//' 
//' @param pnull R function (cdf)
//' @param rnull R function (generate data under null hypothesis)
//' @param qnull R function (quantiles under null hypothesis)
//' @param ralt  R function to generate data under alternative
//' @param param_alt parameters of ralt
//' @param  w (Optional) function to calculate weights, returns -99 if no weights
//' @param phat  function to estimate parameters from the data
//' @param TS function to calculate test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param B  =c(1000, 1000) Number of simulation runs for power and null distribution
//' @param alpha =0.05, type I error of test 
//' @keywords internal
//' @return A matrix of powers
// [[Rcpp::export]]
Rcpp::NumericMatrix power_cont(
        Rcpp::Function pnull, 
        Rcpp::Function rnull,       
        Rcpp::Function qnull, 
        Rcpp::Function ralt, 
        Rcpp::NumericVector param_alt,
        Rcpp::Function w,
        Rcpp::Function phat,  
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,
        Rcpp::IntegerVector B=Rcpp::IntegerVector::create(1000, 1000), 
        const double alpha=0.05) {
  
  int  i, j, k, np=param_alt.size();
  NumericVector x=ralt(param_alt[0]);
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  NumericVector Fx(x.size()), wx(x.size());

  NumericVector TS_data;
  if(typeTS==1) TS_data=TS(x, pnull, phat(x), qnull);
  if(typeTS==2) {
    Rcpp::List res_w = formals_r(Rcpp::_["fun"]=w);
    if(res_w.size()==1) wx=w(x);
    else wx=w(x, phat(x)); 
    TS_data=TS(x, pnull, phat(x), wx);
  }  
  if(typeTS==3) TS_data=TS(x, pnull, phat(x));
  if(typeTS==4) TS_data=TS(x, pnull, phat(x), TSextra);  
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector methods=TS_data.names();
  NumericMatrix tmp(np, nummethods),out(np, nummethods);
  colnames(out) = methods;
  for(i=0;i<B(0);++i) {
     for(j=0;j<np;++j) {
         NumericVector x=ralt(param_alt[j]); 
         TSextra["p"] = phat(x);
         tmp = gof_cont(x, pnull, rnull, qnull, w, phat, TS, typeTS, TSextra, B(1));
         for(k=0;k<nummethods;++k) {
           if(tmp(1,k)<alpha) out(j,k) = out(j,k)+1;
         } 
     }   
  }
  return out/B(0);
}

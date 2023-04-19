#include <Rcpp.h>
#include "TS_disc.h"
#include "nm_calc.h"

using namespace Rcpp;

//' run gof tests for discrete data
//' 
//' @param x an integer vector of counts
//' @param pnull cumulative distribution function under the null hypothesis
//' @param rnull R function (generate data under null hypothesis)
//' @param vals numeric vector of values of discrete random variables.
//' @param phat function to estimate parameters
//' @param TS function that calculates test statistics
//' @param rate =0 rate of Poisson if sample size is random, 0 otherwise.
//' @param B (=5000) Number of simulation runs  
//' @keywords internal
//' @return A matrix of numbers
// [[Rcpp::export]]
NumericMatrix gof_disc(Rcpp::IntegerVector x, 
                       Rcpp::Function pnull, 
                       Rcpp::Function rnull, 
                       Rcpp::NumericVector vals,
                       Rcpp::Function phat, 
                       Rcpp::Function TS,
                       double rate=0,
                       int B=5000) {
  int k=x.size(), n, i, j;
  n = sum(x);
  NumericMatrix nm=nm_calc(n);
  NumericVector Fx(k);
  for(i=0;i<k;++i) Fx(i)=(1.0+i)/(1.0+i);
  NumericVector TS_data=TS(x, Fx, nm, vals); 
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector allMethods=TS_data.names();
  NumericVector TS_sim(nummethods),pvals(nummethods);
  IntegerVector xsim(k);
  NumericMatrix out(2, nummethods);
  colnames(out) = allMethods;
  NumericVector p=phat(x);
  NumericVector pn(vals.size());
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List res = formals_r(Rcpp::_["fun"]=pnull);
  if(res.size()==0) pn=pnull();
  else pn=pnull(p);
  TS_data=TS(x, pn, nm, vals);
/* run simulation to find null distribution */
  for(i=0;i<B;++i) {
    Rcpp::List resr = formals_r(Rcpp::_["fun"]=rnull);
    if(resr.size()==0) xsim=rnull();
    else xsim=rnull(p);
    NumericVector psim(vals.size());
    if(res.size()!=0) {
      psim=phat(xsim); 
      pn=pnull(psim);
    }
    if(rate>0) nm=nm_calc(sum(xsim));
    TS_sim=TS(xsim, pn, nm, vals);
    for(j=0;j<nummethods;++j) {
      if(TS_data(j)<TS_sim(j)) pvals(j)=pvals(j)+1;
    }
  }
/* record test statistics and p values */  
  for(j=0;j<nummethods;++j) {
      out(0, j)=TS_data(j);
      out(1, j)=pvals(j)/B;
  }
 
  return out;

}

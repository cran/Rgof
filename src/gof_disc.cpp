#include <Rcpp.h>
#include "TS_disc.h"
#include "nm_calc.h"

using namespace Rcpp;

//' run gof tests for discrete data
//' 
//' @param x an integer vector of counts
//' @param pnull cumulative distribution function under the null hypothesis
//' @param vals numeric vector of values of discrete random variables.
//' @param rnull R function (generate data under null hypothesis)
//' @param phat function to estimate parameters
//' @param rate =0 rate of Poisson if sample size is random, 0 otherwise.
//' @param B (=5000) Number of simulation runs  
//' @param doMethod List methods to include
//' @keywords internal
//' @return A matrix of numbers
// [[Rcpp::export]]
NumericMatrix gof_disc(Rcpp::IntegerVector x, 
                       Rcpp::Function pnull, 
                       Rcpp::NumericVector vals,                        
                       Rcpp::Function rnull, 
                       Rcpp::Function phat, 
                       double rate=0,
                       int B=5000,
             Rcpp::CharacterVector doMethod=CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1")) {
  int const nummethods=9;
  int k=x.size(), n, i, j;

  NumericVector  TS_data(nummethods), TS_sim(nummethods),pvals(nummethods);
  IntegerVector xsim(k);
  NumericMatrix out(2, nummethods);
  colnames(out) = CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1");
  n = sum(x);
  NumericMatrix nm=nm_calc(n);
  NumericVector p=phat(x);
  NumericVector pn(vals.size());
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List res = formals_r(Rcpp::_["fun"]=pnull);
  if(res.size()==0) pn=pnull();
  else pn=pnull(p);
  TS_data=TS_disc(x, pn, nm, vals, doMethod);
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
    TS_sim=TS_disc(xsim, pn, nm, vals, doMethod);
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

#include <Rcpp.h>
#include "TS_cont.h"

using namespace Rcpp;

//' run gof tests for continuous data
//' 
//' @param x A numeric vector of data
//' @param pnull R function (cdf)
//' @param phat  function to set or estimate parameters of pnull
//' @param rnull R function (generate data under null hypothesis)
//' @param qnull R function (quantiles under null hypothesis)
//' @param B (=5000) Number of simulation runs  
//' @param doMethod List methods to include
//' @keywords internal
//' @return A matrix of numbers
// [[Rcpp::export]]
Rcpp::NumericMatrix gof_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::Function phat,
        Rcpp::Function rnull, 
        Rcpp::Function qnull, 
        int B=5000,
        Rcpp::CharacterVector doMethod=Rcpp::CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1")
        ) {
  int const nummethods=9;
  int n=x.size(), i, j;
  NumericVector xsim(n), TS_data(nummethods),TS_sim(nummethods), pvals(nummethods), Fx(x.size());
  NumericMatrix out(2, nummethods);
  colnames(out) = CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1");
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  NumericVector p=phat(x);
  Rcpp::List res_pnull = formals_r(Rcpp::_["fun"]=pnull);
  if(res_pnull.size()==1) Fx=pnull(x);
  else Fx=pnull(x, p);
  TS_data=TS_cont(x, Fx, p, qnull, doMethod);
  Rcpp::List res_rnull = formals_r(Rcpp::_["fun"]=rnull);
  for(i=0;i<B;++i) {
    if(res_rnull.size()==0) xsim=rnull();
    else xsim=rnull(p);
    NumericVector psim=phat(xsim);
    if(res_pnull.size()==1) Fx=pnull(xsim);
    else Fx=pnull(xsim, psim);
    TS_sim=TS_cont(xsim, Fx, psim, qnull, doMethod);
    for(j=0;j<nummethods;++j) {
      if(TS_data(j)<TS_sim(j)) pvals(j)=pvals(j)+1;
    }
  }
  for(j=0;j<nummethods;++j) {
      out(0, j)=TS_data(j);
      out(1, j)=pvals(j)/B;
  }

  return out;
}

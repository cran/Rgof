#include <Rcpp.h>

using namespace Rcpp;

//' run gof tests for continuous data
//' 
//' @param x A numeric vector of data
//' @param pnull R function (cdf)
//' @param rnull R function (generate data under null hypothesis)
//' @param qnull R function (quantiles under null hypothesis)
//' @param phat  function to set or estimate parameters of pnull 
//' @param TS function that calculates test statistics
//' @param B (=5000) Number of simulation runs 
//' @keywords internal
//' @return A matrix of numbers
// [[Rcpp::export]]
Rcpp::NumericMatrix gof_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::Function rnull, 
        Rcpp::Function qnull,
        Rcpp::Function phat, 
        Rcpp::Function TS,
        int B=5000) {
  int n=x.size(), i, j;
  NumericVector Fx(x.size());
  NumericVector psim=phat(x);
  for(i=0;i<n;++i) Fx(i)=(1.0+i)/(2.0+i);
  NumericVector TS_data=TS(x, Fx, 1.0, qnull);
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector allMethods=TS_data.names();
  NumericVector xsim(n), TS_sim(nummethods), pvals(nummethods);
  NumericMatrix out(2, nummethods);
  colnames(out) = allMethods;
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  NumericVector p=phat(x);
  Rcpp::List res_pnull = formals_r(Rcpp::_["fun"]=pnull);
  if(res_pnull.size()==1) Fx=pnull(x);
  else Fx=pnull(x, p);
  TS_data=TS(x, Fx, p, qnull);
  Rcpp::List res_rnull = formals_r(Rcpp::_["fun"]=rnull);
  for(i=0;i<B;++i) {
    if(res_rnull.size()==0) xsim=rnull();
    else xsim=rnull(p);
    psim=phat(xsim);
    if(res_pnull.size()==1) Fx=pnull(xsim);    
    else Fx=pnull(xsim, psim);    
    TS_sim=TS(xsim, Fx, psim, qnull);
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

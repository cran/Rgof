#include <Rcpp.h>
#include "TS_disc.h"

using namespace Rcpp;

//' run gof tests for discrete data
//' 
//' @param x an integer vector of counts
//' @param pnull cumulative distribution function under the null hypothesis
//' @param rnull R function (generate data under null hypothesis)
//' @param vals numeric vector of values of discrete random variables.
//' @param phat function to estimate parameters
//' @param TS function that calculates test statistics
//' @param typeTS type of test statistic
//' @param TSextra list passed to TS, if desired
//' @param rate =0, rate of sample size, if random
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
                       int typeTS, 
                       Rcpp::List TSextra,
                       double rate=0.0,
                       int B=5000) {
  int k=x.size(), i, j;
  NumericVector TS_data;
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List res_TS = formals_r(Rcpp::_["fun"]=TS);
  if(typeTS<=1) TS_data=TS(x, pnull, phat(x), vals); 
  if(typeTS==2) TS_data=TS(x, pnull, phat(x), vals, TSextra); 
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector allMethods=TS_data.names();
  NumericVector TS_sim(nummethods),pvals(nummethods);
  IntegerVector xsim(k);
  NumericMatrix out(2, nummethods);
  colnames(out) = allMethods;

/* run simulation to find null distribution */
  for(i=0;i<B;++i) {
    Rcpp::List resr = formals_r(Rcpp::_["fun"]=rnull);
    if(resr.size()==0) xsim=rnull();
    else xsim=rnull(phat(x));
    NumericVector psim(vals.size());
    if(resr.size()!=0) psim=phat(xsim); 
    if(typeTS<=1) TS_sim=TS(xsim, pnull, psim, vals); 
    if(typeTS==2) TS_sim=TS(xsim, pnull, psim, vals, TSextra); 
    for(j=0;j<nummethods;++j) {
      if(TS_data(j)<TS_sim(j)) pvals(j)=pvals(j)+1;
    }
  }
/* record test statistics and p.values */  
  for(j=0;j<nummethods;++j) {
      out(0, j)=TS_data(j);
      out(1, j)=pvals(j)/B;
  }
 
  return out;

}

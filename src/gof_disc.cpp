#include <Rcpp.h>
#include "TS_disc.h"
#include "teststatistics.h"
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
  List dta=List::create(Named("x") =x, Named("vals")=vals);
  TS_data=calcTS(dta, pnull, phat(x), TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector allMethods=TS_data.names();
  NumericVector TS_sim(nummethods),pvals(nummethods);
  IntegerVector xsim(k);
  NumericMatrix out(2, nummethods);
  colnames(out) = allMethods;

/* run simulation to find null distribution */
  for(i=0;i<B;++i) {
    NumericVector p=phat(x);
    if(std::abs(p(0)+99)<0.01) xsim=rnull();
    else xsim=rnull(p);
    dta["x"]=xsim;
    NumericVector TS_sim=calcTS(dta, pnull, phat(xsim), TS, typeTS, TSextra);
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

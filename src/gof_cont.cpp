#include <Rcpp.h>
using namespace Rcpp;

//' run gof tests for continuous data
//' 
//' @param x A numeric vector of data
//' @param pnull R function (cdf)
//' @param rnull R function (generate data under null hypothesis)
//' @param qnull R function (quantiles under null hypothesis)
//' @param w function to calculate weights, returns -99 if no weights
//' @param phat  function to set or estimate parameters of pnull 
//' @param TS function that calculates test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param B (=5000) Number of simulation runs 
//' @keywords internal
//' @return A matrix of numbers
// [[Rcpp::export]]
Rcpp::NumericMatrix gof_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::Function rnull, 
        Rcpp::Function qnull,
        Rcpp::Function w,
        Rcpp::Function phat, 
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,
        int B=5000) {

  int n=x.size(), i, j;
  NumericVector Fx(n), wx(n);
  NumericVector p=phat(x), psim=phat(x);
/* Check the number of arguments of various functions and run them accordingly */  
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List res_pnull = formals_r(Rcpp::_["fun"]=pnull);
  if(res_pnull.size()==1) Fx=pnull(x);
  else Fx=pnull(x, p);
  Rcpp::List res_w = formals_r(Rcpp::_["fun"]=w);
  if(res_w.size()==1) wx=w(x);
  else wx=w(x, p);  
  Rcpp::List res_rnull = formals_r(Rcpp::_["fun"]=rnull);
/* Find test statistics for the data */  
  NumericVector TS_data;
  if(typeTS==1) TS_data=TS(x, pnull, p, qnull);
  if(typeTS==2) TS_data=TS(x, pnull, p, wx);  
  if(typeTS==3) TS_data=TS(x, pnull, p);
  if(typeTS==4) TS_data=TS(x, pnull, p, TSextra);
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector methods=TS_data.names();
  NumericVector xsim(n), TS_sim(nummethods), pvals(nummethods);
  NumericMatrix out(2, nummethods);
  colnames(out) = methods;
  for(i=0;i<B;++i) {
    if(res_rnull.size()==0) xsim=rnull();
    else xsim=rnull(p);
    psim=phat(xsim);
    if(typeTS==1) TS_sim=TS(xsim, pnull, psim, qnull);
    if(typeTS==2) {
      if(res_w.size()==1) wx=w(xsim);
      else wx=w(xsim, psim); 
      TS_sim=TS(xsim, pnull, psim, wx);
    }  
    if(typeTS==3) TS_sim=TS(xsim, pnull, psim);
    if(typeTS==4) {
      TSextra["p"]=psim; 
      TS_sim=TS(xsim, pnull, psim, TSextra);    
    }   
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

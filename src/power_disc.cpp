#include <Rcpp.h>
#include "teststatistics.h"
using namespace Rcpp;

//' find power of gof tests for continuous data
//' 
//' @param pnull R function (cdf)
//' @param rnull R function (generate data under null hypothesis)
//' @param vals values of discrete distribution
//' @param ralt  R function to generate data under alternative
//' @param param_alt parameters of ralt
//' @param phat  function to estimate parameters from the data
//' @param TS function to calculate test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param B  =1000 Number of simulation runs
//' @keywords internal
//' @return A matrix of powers
// [[Rcpp::export]]
Rcpp::List power_disc(
        Rcpp::Function pnull, 
        Rcpp::Function rnull,       
        Rcpp::NumericVector vals, 
        Rcpp::Function ralt, 
        Rcpp::NumericVector param_alt,
        Rcpp::Function phat,  
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,
        const int B=1000) {

  int  i, l, m, np=param_alt.size();
  IntegerVector x=ralt(param_alt[0]);
  NumericVector p=phat(x);
  NumericVector TS_data;
  List dta=List::create(Named("x") =x, Named("vals")=vals);
  int withest=0;
  if(std::abs(p(0)+99)>0.001) withest=1;
  
  TS_data=calcTS(dta, pnull, phat(x), TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  NumericMatrix realdata(B, nummethods), simdata(B*np, 1+nummethods);
  Rcpp::CharacterVector methods=TS_data.names();
  int cn=-1;
  for(l=0;l<B;++l) {
     if(withest==0) x=rnull();
     else x=rnull(p);
     p=phat(x);
     dta["x"]=x;
     TS_data = calcTS(dta, pnull, phat(x), TS, typeTS,  TSextra);
     for(i=0;i<nummethods;++i) realdata(l,i)=TS_data(i);
     for(m=0;m<np;++m) {
         ++cn;
         dta["x"]=ralt(param_alt[m]); 
         TS_data = calcTS(dta, pnull, phat(x),TS, typeTS, TSextra);
         simdata(cn,0)=param_alt(m);
         for(i=0;i<nummethods;++i) simdata(cn,i+1)=TS_data(i);
     }   
  }
  return List::create(Named("Data")=realdata, 
                      Named("Sim")=simdata);
}

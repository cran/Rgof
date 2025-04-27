#include <Rcpp.h>
#include "calcTS.h"
using namespace Rcpp;

//' find power of gof tests for continuous data
//' 
//' @param rnull R function (generate data under null hypothesis)
//' @param vals values of discrete random variable
//' @param ralt  R function to generate data under alternative
//' @param param_alt parameters of ralt
//' @param TS function to calculate test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param B  =1000 Number of simulation runs
//' @keywords internal
//' @return A matrix of powers
// [[Rcpp::export]]
Rcpp::List gof_power_C(
        Rcpp::Function rnull,
        Rcpp::NumericVector vals,
        Rcpp::Function ralt, 
        Rcpp::NumericVector param_alt,
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,
        const int B=1000) {

  int  i, l, m, np=param_alt.size();
  List dta=List::create(Named("x"),Named("vals")=vals);
  dta["x"]=ralt(param_alt[0]);
  NumericVector TS_data = calcTS(dta, TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
  Rcpp::List resralt = formals_r(Rcpp::_["fun"]=ralt);
  Rcpp::Function phat=TSextra["phat"];
  NumericVector p=phat(dta["x"]), tmp1;
  IntegerVector tmp2;
  NumericMatrix realdata(B, nummethods), simdata(B*np, 1+nummethods);
  Rcpp::CharacterVector methods=TS_data.names();
  int cn=-1;
  for(l=0;l<B;++l) {
     if(TSextra["Continuous"]) {
       if(resrnull.size()==0) tmp1=rnull();
       else tmp1=rnull(p);
       dta["x"]=tmp1;
     }
     else {
       if(resrnull.size()==0) tmp2=rnull();
       else tmp2=rnull(p);
       dta["x"]=tmp2;
       dta["vals"]=vals;
     }
     TS_data = calcTS(dta, TS, typeTS, TSextra); 
     for(i=0;i<nummethods;++i) realdata(l,i)=TS_data(i);
     for(m=0;m<np;++m) {
         ++cn;
         if(resralt.size()==1) dta["x"]=ralt(param_alt[m]);
         else dta["x"]=ralt(param_alt[m], p);
         TS_data = calcTS(dta, TS, typeTS, TSextra);
         simdata(cn,0)=param_alt(m);
         for(i=0;i<nummethods;++i) simdata(cn,i+1)=TS_data(i);
     }   
  }
  return List::create(Named("Data")=realdata, 
                      Named("Sim")=simdata);
}

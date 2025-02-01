#include <Rcpp.h>
#include "Cpporder.h"
using namespace Rcpp;

//' Find test statistics for continuous data
//' 
//' @param x A numeric vector.
//' @param pnull cdf.
//' @param param parameters for pnull  in case of parameter estimation.
//' @param qnull An R function, the quantile function under the null hypothesis.
//' @return A numeric vector with test statistics
// [[Rcpp::export]]
Rcpp::NumericVector TS_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::NumericVector param, 
        Rcpp::Function qnull) {
  Rcpp::CharacterVector methods=Rcpp::CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1");
  int const nummethods=methods.size();
  int n=x.size(), i;
  NumericVector TS(nummethods), b1(n), b2(n);
  double m, tmp, tmp1;
  TS.names() =  methods;

  /* some setup */
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List respnull = formals_r(Rcpp::_["fun"]=pnull);
  NumericVector Fx(n);
  if(respnull.size()==1) Fx = pnull(x);
  else Fx = pnull(x, param);
  Rcpp::List resqnull = formals_r(Rcpp::_["fun"]=qnull);
  NumericVector qtmp(n);
  if(resqnull.size()==1) qtmp = qnull(0.5);
  else qtmp = qnull(0.5, 0);
  
  /*  sort data */

  x = Cpporder(x, x);
  Fx = Cpporder(Fx, Fx);

  /*  Kolmogorov-Smirnov and Kuiper*/

    TS(0)=0;
    double mx=0.0, Mx=0.0;
    for(i=0;i<n;++i) {
      tmp = Fx[i]-double(i)/n;
      if(tmp<mx) mx=tmp;
      if(tmp>Mx) Mx=tmp;
      tmp = double(i+1)/n-Fx[i];
      if(tmp<mx) mx=tmp;
      if(tmp>Mx) Mx=tmp;
    }
    if(std::abs(mx)>std::abs(Mx)) TS(0)=std::abs(mx);
    else TS(0)=std::abs(Mx);
    TS(1) = std::abs(mx)+std::abs(Mx);

  /* Anderson-Darling */

 
    tmp=0.0;
    for(i=0;i<n;++i) {
      tmp=tmp + (2*i+1)*(log(Fx[i])+log(1.0-Fx[n-i-1]));        
    }
    TS(2)=-n-tmp/n;
  
  /* Cramer-von Mises and/or Watson*/

    tmp=0.0;
    tmp1=0.0;
    for(i=0;i<n;++i) {
      tmp = tmp + ((2*i+1)/2.0/n-Fx[i])*((2*i+1)/2.0/n-Fx[i]);
      tmp1 = tmp1 + Fx[i];
    }
    TS(3) = 1/(12.0*n)+tmp;
    TS(4) = TS(3) - n*(tmp1/n-0.5)*(tmp1/n-0.5);

  /* Zhang's Methods */

    TS(5) = 0.0;
    TS(6) = 0.0;
    TS(7) = 0.0;
    for(i=0; i<n; ++i) {
      if(Fx[i]<1) {
        m = i+0.5;
        TS(5) = TS(5) - log(Fx[i])/(n-m) - log(1-Fx[i])/m;
        tmp =  m*log(m/n/Fx[i]) + (n-m)*log((n-m)/n/(1-Fx[i]));
        if(tmp>TS(6)) TS(6)=tmp;
        tmp = log( (1/Fx[i]-1)/((n-0.5)/(i+0.25)-1) );
        TS(7) = TS(7) + tmp*tmp;
      }  
    }
  
    NumericVector a1(n+1),a2(n);
    for(i=0;i<=n;++i) {
      a1[i]=double(i)/n;
      if(a1[i]<1e-10) a1[i]=1e-10;
      if(a1[i]>1-1e-10) a1[i]=1-1e-10;
    }  
    if(resqnull.size()==1) b1=qnull(a1);
    else b1=qnull(a1, param);
    for(i=0;i<n;++i) a2[i]=(2.0*i+1.0)/2.0/n;
    if(resqnull.size()==1) b2=qnull(a2);
    else b2=qnull(a2, param);
    
    TS(8) = 0.0;
    for(i=0;i<n;++i) {
     TS(8) = TS(8) + std::abs(b1(i)-x[i]);
     TS(8) = TS(8) + 4*std::abs(b2(i)-x[i]);
     TS(8) = TS(8) + std::abs(b1(i+1)-x[i]);
   }
   TS(8) = TS(8)/6.0/n;

  return TS;
} 

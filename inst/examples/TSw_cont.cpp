#include <Rcpp.h>
#include "Cpporder.h"
using namespace Rcpp;

//' Find test statistics for continuous data with weights
//' 
//' @param x A numeric vector.
//' @param Fx numeric vector of cdf probabilities.
//' @param w numeric vector of weights
//' @keywords internal
//' @return A numeric vector with test statistics
// [[Rcpp::export]]
Rcpp::NumericVector TSw_cont(
        Rcpp::NumericVector x, 
        Rcpp::NumericVector Fx,
        Rcpp::NumericVector w) {
  
  Rcpp::CharacterVector methods=Rcpp::CharacterVector::create("KS", "K", "CvM", "AD");
  int const nummethods=methods.size();
  int n=x.size(), i;
  NumericVector TS(nummethods), edf(n);
  double tmp, tmp1;
  TS.names() =  methods;

  /*  sort data */

  w = Cpporder(w, x);
  Fx = Cpporder(Fx, x);
  x = Cpporder(x, x);

  /* weighted edf */

    edf(0) = w(0);
    for(i=1;i<n;++i) edf(i)=edf(i-1)+w(i);
    double sw=edf(n-1);
    for(i=0;i<n;++i) edf(i)=edf(i)/sw;

  /*  Kolmogorov-Smirnov and Kuiper*/

    TS(0)=0;
    double mx=0.0, Mx=0.0;
    for(i=0;i<n;++i) {
      tmp = edf(i)-Fx(i);
      if(std::abs(tmp)>std::abs(mx)) mx=std::abs(tmp);
      if(i>1) {
        tmp = Fx(i)-edf(i-1);
        if(std::abs(tmp)>std::abs(Mx)) Mx=std::abs(tmp);
      }  
    }
    if(std::abs(mx)>std::abs(Mx)) TS(0)=std::abs(mx);
    else TS(0)=std::abs(Mx);
    TS(1) = std::abs(mx)+std::abs(Mx);

    /* Cramer-von Mises*/

    TS(2) = n/3.0;  
    for(i=0;i<n;++i) {
      if(i==0) tmp = edf(0)/2.0;
      else tmp = (edf(i)+edf(i-1))/2.0;
      TS(2) = TS(2) + 
        n*w(i)/sw*((Fx(i)-tmp)*(Fx(i)-tmp) - tmp*tmp);
    }

  /* Anderson-Darling */
 
    TS(3)=0.0;
    for(i=0;i<n;++i) {
      if(i==0) {
        tmp = -edf(0)*edf(0);  
        tmp1 = (edf(0)-1)*(edf(0)-1)-1;
      }  
      else {
        tmp = edf(i-1)*edf(i-1)-edf(i)*edf(i);  
        tmp1 = (edf(i)-1)*(edf(i)-1)-(edf(i-1)-1)*(edf(i-1)-1);        
      }
      if(std::abs(Fx[i]-0.5)<0.5) 
         TS(3) = TS(3) + tmp*log(Fx[i]) + tmp1*log(1.0-Fx[i]); 
    }
    TS(3)=-n+n*TS(3);
  
  return TS;
} 

#include <Rcpp.h>
using namespace Rcpp;

//' Find test statistics for continuous data
//' 
//' @param x A numeric vector.
//' @param Fx numeric vector of cdf probabilities.
//' @param param parameters for pnull
//' @param qnull An R function, the quantile function under the null hypothesis.
//' @param doMethod A character vector of methods to include
//' @keywords internal
//' @return A numeric vector with test statistics
// [[Rcpp::export]]
Rcpp::NumericVector TS_cont(
        Rcpp::NumericVector x, 
        Rcpp::NumericVector Fx,
        Rcpp::NumericVector param, 
        Rcpp::Function qnull,
        Rcpp::CharacterVector doMethod=CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1")) {

  int const nummethods=9;
  int n=x.size(), i;
  NumericVector TS(nummethods), b1(n), b2(n);
  double m, tmp, tmp1;

  TS.names() =  CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1");

  /*  sort data */

  std::sort(x.begin(), x.end());
  std::sort(Fx.begin(), Fx.end());

  /*  Kolmogorov-Smirnov and Kuiper*/

  LogicalVector H=in(CharacterVector::create("KS", "K"), doMethod);
  if( (H[0]==TRUE) || (H[1]==TRUE) ) {
    TS(0)=0;
    double mx=0.0, Mx=0.0;
    for(i=0;i<n;++i) {
      tmp = Fx[i]-double(i)/n;
      if(tmp<0 && std::abs(tmp)>std::abs(mx)) mx=std::abs(tmp);
      if(tmp>0 && std::abs(tmp)>std::abs(Mx)) Mx=std::abs(tmp);
      tmp = double(i+1)/n-Fx[i];
      if(tmp<0 && std::abs(tmp)>std::abs(mx)) mx=std::abs(tmp);
      if(tmp>0 && std::abs(tmp)>std::abs(Mx)) Mx=std::abs(tmp);      
    }
    if(std::abs(mx)>std::abs(Mx)) TS(0)=std::abs(mx);
    else TS(0)=std::abs(Mx);
    if(H[1]==TRUE) TS(1)=std::abs(Mx)+std::abs(mx);
    if(in(CharacterVector::create("KS"), doMethod)[0]==FALSE) TS(0)=0.0;
  }

  /* Anderson-Darling */

  if(in(CharacterVector::create("AD"), doMethod)[0]==TRUE) {
    tmp=0.0;
    for(i=0;i<n;++i) {
      tmp=tmp + (2*i+1)*(log(Fx[i])+log(1.0-Fx[n-i-1]));        
    }
    TS(2)=-n-tmp/n;
    if(in(CharacterVector::create("AD"), doMethod)[0]==FALSE) TS(2)=0.0;
  }

  /* Cramer-von Mises and/or Wilson*/

  H=in(CharacterVector::create("CvM", "W"), doMethod);
  if( (H[0]==TRUE) || (H[1]==TRUE) ) {
    tmp=0.0;
    tmp1=0.0;
    for(i=0;i<n;++i) {
      tmp = tmp + ((2*i+1)/2.0/n-Fx[i])*((2*i+1)/2.0/n-Fx[i]);
      tmp1 = tmp1 + Fx[i];
    }
    TS(3) = 1/(12.0*n)+tmp;
    TS(4) = TS(3) - n*(tmp1/n-0.5)*(tmp1/n-0.5);
    
    if(in(CharacterVector::create("CvM"), doMethod)[0]==FALSE) TS(3)=0.0;
    if(in(CharacterVector::create("W"), doMethod)[0]==FALSE) TS(4)=0.0;
  }

  /* Zhang's Methods */

  H=in(CharacterVector::create("ZA", "ZK", "ZC"), doMethod);
  if( (H[0]==TRUE) || (H[1]==TRUE) || (H[2]==TRUE) ) {
    TS(5) = 0.0;
    TS(6) = 0.0;
    TS(7) = 0.0;
    for(i=0; i<n; ++i) {
      m = i+0.5;
      TS(5) = TS(5) - log(Fx[i])/(n-m) - log(1-Fx[i])/m;
      tmp =  m*log(m/n/Fx[i]) + (n-m)*log((n-m)/n/(1-Fx[i]));
      if(tmp>TS(6)) TS(6)=tmp;
      tmp = log( (1/Fx[i]-1)/((n-0.5)/(i+0.25)-1) );
      TS(7) = TS(7) + tmp*tmp;
    }
    if(in(CharacterVector::create("ZA"), doMethod)[0]==FALSE) TS(5)=0.0;
    if(in(CharacterVector::create("ZK"), doMethod)[0]==FALSE) TS(6)=0.0;
    if(in(CharacterVector::create("ZC"), doMethod)[0]==FALSE) TS(7)=0.0;
  }

  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List res = formals_r(Rcpp::_["fun"]=qnull);
  NumericVector qtmp(n);
  if(res.size()==1) qtmp = qnull(0.5);
  else qtmp = qnull(0.5, 0);
  if( (in(CharacterVector::create("Wassp1"), doMethod)[0]==TRUE) && (qtmp[0]!=-99)   ) {
    
    NumericVector a1(n+1),a2(n);
    for(i=0;i<=n;++i) {
      a1[i]=double(i)/n;
      if(a1[i]<1e-10) a1[i]=1e-10;
      if(a1[i]>1-1e-10) a1[i]=1-1e-10;
    }  
    if(res.size()==1) b1=qnull(a1);
    else b1=qnull(a1, param);
    for(i=0;i<n;++i) a2[i]=(2.0*i+1.0)/2.0/n;
    if(res.size()==1) b2=qnull(a2);
    else b2=qnull(a2, param);
    
    TS(8) = 0.0;
    for(i=0;i<n;++i) {
     TS(8) = TS(8) + std::abs(b1(i)-x[i]);
     TS(8) = TS(8) + 4*std::abs(b2(i)-x[i]);
     TS(8) = TS(8) + std::abs(b1(i+1)-x[i]);
   }
   TS(8) = TS(8)/6.0/n;
  
  }

  return TS;
} 

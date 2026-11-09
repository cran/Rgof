#include <Rcpp.h>
using namespace Rcpp;

//' Find test statistics for discrete data
//' 
//' @param x An integer vector.
//' @param pnull cdf.
//' @param param parameters for pnull in case of parameter estimation.
//' @param vals A numeric vector with the values of the discrete rv.
//' @param w weights 
//' @keywords internal
//' @return A vector with test statistics
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector TSw_disc(
        Rcpp::IntegerVector x, 
        Rcpp::Function pnull,
        Rcpp::NumericVector param,
        Rcpp::NumericVector vals,  
        Rcpp::NumericVector w) {
    
  Rcpp::CharacterVector methods=CharacterVector::create("KS", "K", "CvM", "AD");    
  int const nummethods=methods.size();
  int k=x.size(),  i;
  NumericVector TS(nummethods), edf(k);
  double n, tmp;
  TS.names() =  methods;

  /* some setup */
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List respnull = formals_r(Rcpp::_["fun"]=pnull);
  NumericVector Fx(k);
  if(respnull.size()==0) Fx = pnull();
  else Fx = pnull(param);
  
  /*  Find sample size, cumulative sum of x and a vector used in various calculations*/
  
  n=0.0;
  for(i=0;i<k;++i) n = n + x[i]*w[i];

  /* Edf */
  
  edf[0] = x[0]*w[0]/n;
  for(i=1;i<k;++i) edf[i] = edf[i-1] + x[i]*w[i]/n;
  
  /*  Kolmogorov-Smirnov and Kuiper*/
    double mx = 0;
    double Mx = 0;
    for(i=0;i<k-1;++i) {
      tmp = edf[i] - Fx[i];
      if(tmp<0 && std::abs(tmp)>std::abs(mx)) mx=std::abs(tmp);
      if(tmp>0 && std::abs(tmp)>std::abs(Mx)) Mx=std::abs(tmp);      
    }
    if(std::abs(mx)>std::abs(Mx)) TS(0)=std::abs(mx);
    else TS(0)=std::abs(Mx);
    TS(1)=Mx+mx; 
 
 /* Cramer-von Mises and  Anderson-Darling, with the weights
    t(j)=(p(j)+p(j+1))/2 of Choulakian, Lockhart and Stephens (1994) */
 
    NumericVector p(k), t(k);
    for(i=0;i<k;++i) {
      if(i==0) p[0]=Fx[0];
      else p[i]=Fx[i]-Fx[i-1];
    }
    for(i=0;i<k;++i) {
      if(i<k-1) t[i]=(p[i]+p[i+1])/2.0;
      else t[i]=(p[k-1]+p[0])/2.0;
    }
    TS(2) = 0.0;
    TS(3) = 0.0;
    for(i=0;i<k;++i) {
      tmp = (edf[i]-Fx[i])*(edf[i]-Fx[i]);
      TS(2) = TS(2) + tmp*t[i];
      if( (i<k-1)&&(Fx[i]>0)&&(Fx[i]<1) ) 
        TS(3) = TS(3) + tmp/Fx[i]/(1-Fx[i])*t[i]; 
    }
    TS(2) = n*TS(2);
    TS(3) = n*TS(3);
 
    return TS;
}

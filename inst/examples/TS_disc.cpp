#include <Rcpp.h>
using namespace Rcpp;

//' Find test statistics for discrete data
//' 
//' @param x An integer vector.
//' @param Fx A numeric vector of cumulative probabilities.
//' @param vals A numeric vector with the values of the discrete rv.
//' @keywords internal
//' @return A vector with test statistics
// [[Rcpp::export]]
NumericVector TS_disc(IntegerVector x, 
                      NumericVector Fx,  
                      NumericVector vals) {
    
  Rcpp::CharacterVector methods=CharacterVector::create("KS", "K", "AD", "CvM", "W", "Wassp1");    
  int const nummethods=methods.size();
  int k=x.size(), n, i;
  NumericVector TS(nummethods), ecdf(k);
  NumericMatrix logF(k, 4);
  IntegerVector cumx(k), cumx1(k);
  double tmp, tmp1;
  TS.names() =  methods;

  /*  Find sample size, cumulative sum of x and a vector used in various calculations*/
  
  n=0;
  cumx[0] = x[0];
  cumx1[0] = x[k-1];  
  for(i=0;i<k;++i) {
    n = n + x[i];
    if(i>0) {
      cumx[i] = cumx[i-1] + x[i];
      cumx1[i] = cumx1[i-1] + x[k-1-i];
    }  
  }

  /* Ecdf, distribution function and its logs evaluated at data vals*/
  
  double mFx=0.0;
  for(i=0;i<k;++i) if(Fx(i)>mFx && Fx(i)<1) mFx=Fx(i);  
  
  for(i=0;i<k;++i) {
    logF(i, 0) = log(Fx[i]);
    if(Fx[i]<1) {
       logF(i, 1) = log(1.0-Fx[i]);
    }   
    else {
      tmp = (9.0+mFx)/10.0;
      logF(i, 1) = log(1.0-tmp);
    }  
    if(Fx(k-i-1)<1) logF(i, 3) = log(1-Fx(k-i-1));
    else logF(i, 3) = 0.0;      
  }
   
  ecdf[0] = double(x[0])/double(n);
  for(i=1;i<k;++i) ecdf[i] = ecdf[i-1] + x[i]/double(n);
  
  /*  Kolmogorov-Smirnov and Kuiper*/
    double mx = 0;
    double Mx = 0;
    for(i=0;i<k-1;++i) {
      tmp = ecdf[i]- Fx[i];
      if(tmp<0 && std::abs(tmp)>std::abs(mx)) mx=std::abs(tmp);
      if(tmp>0 && std::abs(tmp)>std::abs(Mx)) Mx=std::abs(tmp);      
      if(i==0) tmp=Fx[i];
      else tmp=Fx[i]-ecdf[i-1];
      if(tmp<0 && std::abs(tmp)>std::abs(mx)) mx=std::abs(tmp);
      if(tmp>0 && std::abs(tmp)>std::abs(Mx)) Mx=std::abs(tmp);      
    }
    if(std::abs(mx)>std::abs(Mx)) TS(0)=std::abs(mx);
    else TS(0)=std::abs(Mx);
    TS(1)=Mx+mx; 
  
  /* Anderson-Darling */
   
    if(Fx[0]<1) tmp = (ecdf[0]-Fx[0])*(ecdf[0]-Fx[0])/(1-Fx[0]);
    else tmp = 0.0;
    for(i=1;i<k-1;++i) {
      if( (Fx[i]>0)&&(Fx[i]<1) )
          tmp = tmp + (ecdf[i]-Fx[i])*(ecdf[i]-Fx[i])/Fx[i]/(1-Fx[i])*(Fx[i]-Fx[i-1]); 
    }
    TS(2)=n*tmp;

  /* Cramer-von Mises */
  
    tmp = (ecdf[0]-Fx[0])*(ecdf[0]-Fx[0])*Fx[0];
    for(i=1;i<k;++i) {
      tmp = tmp + (ecdf[i]-Fx[i])*(ecdf[i]-Fx[i])*(Fx[i]-Fx[i-1]);
    }
    TS(3) = 1/(12.0*n)+n*tmp;

  /* Wilson*/
  
    tmp = double(n)*(4.0*n*n-1)/3.0-4.0*n*cumx[0]*(cumx[0]-1)*Fx[0]-4.0*n*x[0]*Fx[0]+4.0*n*n*x[0]*Fx[0]*Fx[0];
    tmp1 = x[0]*Fx[0];
    for(i=1;i<k;++i) {
      tmp = tmp - 4.0*n*(cumx[i]*(cumx[i]-1)-cumx[i-1]*(cumx[i-1]-1))*Fx[i]-4.0*n*x[i]*Fx[i]+4.0*n*n*x[i]*Fx[i]*Fx[i];
      tmp1 = tmp1 + x[i]*Fx[i];
    }
    TS(4) = 1/(12.0*n)+tmp/4/n/n - n*(tmp1/n-0.5)*(tmp1/n-0.5);
    
    /*Wasserstein */   
    TS(5)=0.0;
    for(i=0;i<k-1;++i) 
      TS(5) = TS(5) + std::abs(cumx[i]/double(n)-Fx[i])*(vals[i+1]-vals[i]);
    

  return TS;
}

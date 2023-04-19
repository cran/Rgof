#include <Rcpp.h>
using namespace Rcpp;

//' Find test statistics for discrete data
//' 
//' @param x An integer vector.
//' @param p A numeric vector of probabilities.
//' @param nm A matrix of pre-calculated (with nm_calc.cpp) numbers needed for Zhangs tests.
//' @param vals A numeric vector with the values of the discrete rv.
//' @keywords internal
//' @return A vector with test statistics
// [[Rcpp::export]]
NumericVector TS_disc(IntegerVector x, 
                      NumericVector p,
                      NumericMatrix nm,  
                      NumericVector vals) {
    
  Rcpp::CharacterVector methods=CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZC", "Wassp1");    
  int const nummethods=methods.size();
  int k=x.size(), n, i;
  NumericVector TS(nummethods), ecdf(k), Fx(k);
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

/* find matrix nm if not provided */

  /* Ecdf, distribution function and its logs evaluated at data vals*/
  
  double mFx=0.0;
  for(i=0;i<k;++i) {
     Fx(i)=p(i);
     if(p(i)>mFx && p(i)<1) mFx=p(i);  
  }   
  for(i=0;i<k;++i) {
    logF(i, 0) = log(Fx[i]);
    if(Fx[i]<1) {
       logF(i, 1) = log(1.0-Fx[i]);
       logF(i, 2) = log(1.0/Fx[i]-1.0);
    }   
    else {
      tmp = (9.0+mFx)/10.0;
      logF(i, 1) = log(1.0-tmp);
      logF(i, 2) = log(1.0/tmp-1.0); 
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
      tmp = ecdf[i] - Fx[i];
      if(tmp<0 && std::abs(tmp)>std::abs(mx)) mx=std::abs(tmp);
      if(tmp>0 && std::abs(tmp)>std::abs(Mx)) Mx=std::abs(tmp);      
    }
    if(std::abs(mx)>std::abs(Mx)) TS(0)=std::abs(mx);
    else TS(0)=std::abs(Mx);
    TS(1)=Mx+mx; 
  
  /* Anderson-Darling */
 
   tmp = (ecdf[0]-Fx[0])*(ecdf[0]-Fx[0])/(1-Fx[0]);
    for(i=1;i<k-1;++i) {
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

  /* Zhang's Methods */
  
    TS(5) = 0.0;
    TS(6) = 0.0;
    tmp=0.0;
    for(i=0; i<k; ++i) {
        if(x[i]>0) {
            TS(5) = TS(5) - (nm(cumx[i], 0)-nm(cumx[i-1], 0))*logF(i, 0) - (nm(cumx[i], 1)-nm(cumx[i-1], 1))*logF(i, 1);
            TS(6) = TS(6) + x[i]*logF(i, 2)*logF(i, 2) - 2*(nm(cumx[i], 2)-nm(cumx[i-1], 2))*logF(i, 2) + (nm(cumx[i], 3)-nm(cumx[i-1], 3));          
        }
    }


    TS(7)=0.0;
    for(i=0;i<k-1;++i) 
      TS(7) = TS(7) + std::abs(cumx[i]/double(n)-Fx[i])*(vals[i+1]-vals[i]);

  return TS;
}

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
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector TS_cont(
        Rcpp::NumericVector x, 
        Rcpp::Function pnull,
        Rcpp::NumericVector param, 
        Rcpp::Function qnull) {
  Rcpp::CharacterVector methods=Rcpp::CharacterVector::create("KS", "K", "AD", "CvM", "W", "ZA", "ZK", "ZC", "Wassp1");
  int const nummethods=methods.size();
  int n=x.size(), i;
  NumericVector TS(nummethods);
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
  
  /*  sort data */

  x = Cpporder(x, x);
  Fx = Cpporder(Fx, Fx);

  /*  Kolmogorov-Smirnov and Kuiper:
      D+ = max(i/n-F(x_(i))), D- = max(F(x_(i))-(i-1)/n), KS = max(D+, D-), 
      Kuiper V = D+ + D- */

    double Dp=0.0, Dm=0.0;
    for(i=0;i<n;++i) {
      tmp = double(i+1)/n-Fx[i];
      if(tmp>Dp) Dp=tmp;
      tmp = Fx[i]-double(i)/n;
      if(tmp>Dm) Dm=tmp;
    }
    if(Dp>Dm) TS(0)=Dp;
    else TS(0)=Dm;
    TS(1) = Dp+Dm;

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
 
 /* Wasserstein: W1 = int_0^1 |Fn^{-1}(u)-F^{-1}(u)| du. 
    On ((i-1)/n, i/n] Fn^{-1}(u)=x_(i). Each such interval is split at 
    u=F(x_(i)), so x_(i)-F^{-1}(u) has constant sign on each piece and 
    int |x_(i)-F^{-1}(u)| du = |x_(i)(b-a) - int_a^b F^{-1}(u) du|.
    The integrals of F^{-1} are found with 3-point Gauss-Legendre quadrature, 
    20 points on the two pieces at the ends, where F^{-1} may be unbounded. */
  NumericVector tmp2;
  if(resqnull.size()==1) tmp2=qnull(0.5);
  else tmp2=qnull(0.5, param);
    
  if(std::abs(tmp2[0]+99)>0.01) {
    const double g3x[3]={-0.7745966692414834, 0.0, 0.7745966692414834};
    const double g3w[3]={0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
    const double g20x[10]={0.0765265211334973, 0.2277858511416451, 0.3737060887154195,
       0.5108670019508271, 0.6360536807265150, 0.7463319064601508, 0.8391169718222188,
       0.9122344282513259, 0.9639719272779138, 0.9931285991850949};
    const double g20w[10]={0.1527533871307258, 0.1491729864726037, 0.1420961093183820,
       0.1316886384491766, 0.1181945319615184, 0.1019301198172404, 0.0832767415767048,
       0.0626720483341091, 0.0406014298003869, 0.0176140071391521};
    /* breakpoints: i/n and F(x_(i)), sorted */
    std::vector<double> t;
    for(i=0;i<=n;++i) t.push_back(double(i)/n);
    for(i=0;i<n;++i) if(Fx[i]>0.0 && Fx[i]<1.0) t.push_back(Fx[i]);
    std::sort(t.begin(), t.end());
    t.erase(std::unique(t.begin(), t.end()), t.end());
    int np=t.size()-1, j, l;
    /* all quadrature nodes, evaluated with one call of qnull */
    std::vector<double> nodes;
    for(j=0;j<np;++j) {
      double a=t[j], b=t[j+1], c=(a+b)/2, h=(b-a)/2;
      if(j==0 || j==np-1) {
        for(l=0;l<10;++l) {nodes.push_back(c-h*g20x[l]); nodes.push_back(c+h*g20x[l]);}
      }
      else for(l=0;l<3;++l) nodes.push_back(c+h*g3x[l]);
    }
    NumericVector u(nodes.begin(), nodes.end()), qu;
    if(resqnull.size()==1) qu=qnull(u);
    else qu=qnull(u, param);
    TS(8) = 0.0;
    int k=0;
    for(j=0;j<np;++j) {
      double a=t[j], b=t[j+1], h=(b-a)/2, Q=0.0;
      if(j==0 || j==np-1) {
        for(l=0;l<10;++l) {Q+=g20w[l]*(qu[k]+qu[k+1]); k+=2;}
      }
      else for(l=0;l<3;++l) {Q+=g3w[l]*qu[k]; ++k;}
      Q=h*Q;
      int ix=int(floor((a+b)/2*n));
      if(ix>n-1) ix=n-1;
      if(ix<0) ix=0;
      TS(8) = TS(8) + std::abs(x[ix]*(b-a)-Q);
    }
  }
  return TS;
} 

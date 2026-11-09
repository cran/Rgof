#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

/* Data-driven Neyman smooth tests and Berk-Jones statistic for discrete data.
   x: counts, Fx: null cdf at the values, k values. 
   Smooth tests: with p(j)=F(v(j))-F(v(j-1)) and the mid-cdf u(j)=F(v(j-1))+p(j)/2,
   the Legendre polynomials in 2u-1 are made orthonormal w.r.t. p (Gram-Schmidt, on 
   the values with p(j)>0), g(1),..,g(K), K=min(10, number of values with p>0 - 1). 
   With b(l)=sum_j x(j)g(l)(v(j))/sqrt(n) and S(k)=b(1)^2+..+b(k)^2,
   NS-B = S(k*) with k* maximizing S(k)-k log(n) (Schwarz rule, Ledwina 1994), 
   NS-A = S(k*) with k* maximizing S(k)-2k (AIC).
   Berk-Jones: BJ = n max_j KL(Fn(v(j)), F(v(j))), j<k with 0<F(v(j))<1, 
   KL(a,b)=a log(a/b)+(1-a)log((1-a)/(1-b)). 
   out[0]=NS-B, out[1]=NS-A, out[2]=BJ */
static void smoothBJ_disc(const std::vector<double>& x, const std::vector<double>& Fx, 
                   double* out) {
  int k=x.size(), i, j, l, q;
  double n=0.0;
  for(i=0;i<k;++i) n += x[i];
  std::vector<double> p(k), u(k);
  for(i=0;i<k;++i) {
    p[i] = (i==0) ? Fx[0] : Fx[i]-Fx[i-1];
    if(p[i]<0) p[i] = 0.0;
    u[i] = ((i==0) ? 0.0 : Fx[i-1]) + p[i]/2.0;
  }
  std::vector<int> J;
  for(i=0;i<k;++i) if(p[i]>0) J.push_back(i);
  int m=J.size(), K=std::min(10, m-1);
  out[0] = 0.0; out[1] = 0.0;
  if(K>=1) {
    std::vector< std::vector<double> > g(K+1, std::vector<double>(m));
    for(j=0;j<m;++j) {
      double t=2.0*u[J[j]]-1.0;
      g[0][j] = 1.0;
      g[1][j] = t;
      for(l=2;l<=K;++l) g[l][j] = ((2*l-1)*t*g[l-1][j]-(l-1)*g[l-2][j])/l;
    }
    /* Gram-Schmidt with inner product sum_j p(j) f(j) g(j), done twice */
    std::vector<bool> ok(K+1, true);
    for(l=0;l<=K;++l) {
      for(int rep=0;rep<2;++rep) {
        for(q=0;q<l;++q) {
          if(!ok[q]) continue;
          double c=0.0;
          for(j=0;j<m;++j) c += p[J[j]]*g[l][j]*g[q][j];
          for(j=0;j<m;++j) g[l][j] -= c*g[q][j];
        }
      }
      double nrm=0.0;
      for(j=0;j<m;++j) nrm += p[J[j]]*g[l][j]*g[l][j];
      nrm = std::sqrt(nrm);
      if(nrm<1e-10) {ok[l]=false; for(j=0;j<m;++j) g[l][j]=0.0;}
      else for(j=0;j<m;++j) g[l][j] /= nrm;
    }
    double S=0.0, bestB=-1e300, bestA=-1e300;
    for(l=1;l<=K;++l) {
      double b=0.0;
      for(j=0;j<m;++j) b += x[J[j]]*g[l][j];
      b = b/std::sqrt(n);
      S += b*b;
      if(S-l*std::log(n)>bestB) {bestB=S-l*std::log(n); out[0]=S;}
      if(S-2.0*l>bestA) {bestA=S-2.0*l; out[1]=S;}
    }
  }
  /* Berk-Jones */
  double cum=0.0, kl, a, f, bj=0.0;
  for(i=0;i<k-1;++i) {
    cum += x[i];
    f = Fx[i];
    if(f<=0.0 || f>=1.0) continue;
    a = cum/n;
    if(a<=0.0) kl = std::log(1.0/(1.0-f));
    else if(a>=1.0) kl = std::log(1.0/f);
    else kl = a*std::log(a/f)+(1.0-a)*std::log((1.0-a)/(1.0-f));
    if(kl>bj) bj=kl;
  }
  out[2] = n*bj;
}


//' Find test statistics for discrete data
//'
//' With Z(j) = Fn(v(j))-F(v(j)), p(j) = F(v(j))-F(v(j-1)) and
//' t(j) = (p(j)+p(j+1))/2, with p(k+1)=p(1), the statistics are
//' KS = max|Z(j)|, Kuiper = max(Z(j),0) + max(-Z(j),0),
//' CvM = n sum Z(j)^2 t(j), Watson = n sum (Z(j)-Zbar)^2 t(j) with
//' Zbar = sum Z(j) t(j), AD = n sum Z(j)^2 t(j)/(F(v(j))(1-F(v(j)))), j<k,
//' as in Choulakian, Lockhart and Stephens (1994), and
//' Wasserstein = sum |Z(j)| (v(j+1)-v(j)).
//' NS-B and NS-A are data-driven Neyman smooth tests, with Legendre polynomials in the
//' mid-cdf made orthonormal w.r.t. the null probabilities and the order chosen by the
//' Schwarz rule (Ledwina 1994) or by AIC, BJ is the Berk-Jones (1979) statistic.
//'
//' @param x An integer vector.
//' @param pnull cdf.
//' @param param parameters for pnull in case of parameter estimation.
//' @param vals A numeric vector with the values of the discrete rv.
//' @return A vector with test statistics
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector TS_disc(
         Rcpp::IntegerVector x,
         Rcpp::Function pnull,
         Rcpp::NumericVector param,
         Rcpp::NumericVector vals) {

  Rcpp::CharacterVector methods=CharacterVector::create("KS", "K", "AD", "CvM", "W", "Wassp1", "NS-B", "NS-A", "BJ");
  int const nummethods=methods.size();
  int k=x.size(), n, i;
  NumericVector TS(nummethods), ecdf(k), Z(k), p(k), t(k);
  double tmp;
  TS.names() =  methods;

  /* some setup */
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List respnull = formals_r(Rcpp::_["fun"]=pnull);
  NumericVector Fx(k);
  if(respnull.size()==0) Fx = pnull();
  else Fx = pnull(param);

  /*  sample size, edf, Z and probabilities */

  n=0;
  for(i=0;i<k;++i) n = n + x[i];
  for(i=0;i<k;++i) {
    if(i==0) {
      ecdf[0] = double(x[0])/double(n);
      p[0] = Fx[0];
    }
    else {
      ecdf[i] = ecdf[i-1] + x[i]/double(n);
      p[i] = Fx[i]-Fx[i-1];
    }
    Z[i] = ecdf[i]-Fx[i];
  }
  /* weights of Choulakian, Lockhart and Stephens (1994) */
  for(i=0;i<k;++i) {
    if(i<k-1) t[i] = (p[i]+p[i+1])/2.0;
    else t[i] = (p[k-1]+p[0])/2.0;
  }

  /*  Kolmogorov-Smirnov and Kuiper. Both Fn and F are step functions
      with jumps at vals, so the extrema are attained at vals */

    double Dp = 0.0, Dm = 0.0;
    for(i=0;i<k;++i) {
      if(Z[i]>Dp) Dp=Z[i];
      if(-Z[i]>Dm) Dm=-Z[i];
    }
    if(Dp>Dm) TS(0)=Dp;
    else TS(0)=Dm;
    TS(1)=Dp+Dm;

  /* Anderson-Darling, last value excluded */
   
    tmp = 0.0;
    for(i=0;i<k-1;++i) {
      if( (Fx[i]>0)&&(Fx[i]<1) )
          tmp = tmp + Z[i]*Z[i]/Fx[i]/(1-Fx[i])*t[i]; 
    }
    TS(2)=n*tmp;

  /* Cramer-von Mises and Watson */
  
    double Zbar = 0.0;
    for(i=0;i<k;++i) Zbar = Zbar + Z[i]*t[i];
    tmp = 0.0;
    double tmp1 = 0.0;
    for(i=0;i<k;++i) {
      tmp = tmp + Z[i]*Z[i]*t[i];
      tmp1 = tmp1 + (Z[i]-Zbar)*(Z[i]-Zbar)*t[i];
    }
    TS(3) = n*tmp;
    TS(4) = n*tmp1;
    
  /*Wasserstein */

    TS(5)=0.0;
    for(i=0;i<k-1;++i)
      TS(5) = TS(5) + std::abs(Z[i])*(vals[i+1]-vals[i]);

  /* Smooth tests and Berk-Jones */

    std::vector<double> xv(k), Fv(k);
    for(i=0;i<k;++i) {xv[i]=x[i]; Fv[i]=Fx[i];}
    double sb[3];
    smoothBJ_disc(xv, Fv, sb);
    TS(6)=sb[0];
    TS(7)=sb[1];
    TS(8)=sb[2];

  return TS;
}

#include <Rcpp.h>
#include <string>
using namespace Rcpp;

//' chi square statistic for discrete data
//' 
//' @param  param parameter for pnull 
//' @param  x  Vector of observations
//' @param  pnull cumulative distribution function under the null hypothesis
//' @param  bins vector of indices of bins
//' @param  formula type of chi square formula to use
//' @param  rate rate of Poisson if sample size is random, 0 otherwise
//' @keywords internal
//' @return chi square statistic
// [[Rcpp::export]]
double chi_stat_disc(
      Rcpp::NumericVector param,    
      Rcpp::IntegerVector x, 
      Rcpp::Function pnull, 
      Rcpp::IntegerVector bins,
      std::string formula,
      double rate=0.0) {
   int n=0, k=bins.size(), i, j, m;
   double tmpO, tmpE;
   NumericVector E(x.size()), p(x.size());  
   n = sum(x);
   if(rate==0) rate=double(n);
   Rcpp::Environment base("package:base");
   Rcpp::Function formals_r = base["formals"];
   Rcpp::List res = formals_r(Rcpp::_["fun"]=pnull);
   if(res.size()==0) p=pnull();
   else p=pnull(param);
   E[0] = rate*p[0];
   for(i=1;i<x.size();++i) E[i] = rate*(p[i]-p[i-1]);
   E = E/max(p);
   double chi = 0.0;
   for(i=0;i<k;++i) {
     tmpO = 0.0;
     tmpE = 0.0;
     if(i<k-2) m = bins[i+1];
     else m = x.size();
     for(j=bins[i];j<m;++j) {
        tmpO = tmpO + x[j];
        tmpE = tmpE + E[j];
     }
     if(formula=="Pearson") chi = chi + (tmpO-tmpE)*(tmpO-tmpE)/tmpE;
     else 
       if(tmpO>0) chi = chi + 2.0*(tmpE-tmpO+tmpO*log(tmpO/tmpE));
     if(m-1==bins[k-1]) i=x.size(); /* skip last bin */   
   }
   return chi;
}

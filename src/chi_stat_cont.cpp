#include <Rcpp.h>
#include "bincounter_cpp.h"
#include <string>
using namespace std;
using namespace Rcpp;

//' chi square statistic for continuous data
//' 
//' @param  param  Vector of parameters of pnull
//' @param  x  Vector of observations
//' @param  pnull Function to calculate probabilities (pmf)
//' @param  bins bins to use
//' @param  formula type of chi square formula to use
//' @param  rate rate of Poisson if sample size is random, 0 otherwise
//' @keywords internal
//' @return chi square statistic
// [[Rcpp::export]]
double chi_stat_cont(
     Rcpp::NumericVector param,
     Rcpp::NumericVector x, 
     Rcpp::Function pnull, 
     Rcpp::NumericVector bins,
     std::string formula,
     double rate=0.0) {
   int n=x.size(), k=bins.size()-1, i;
   NumericVector E(k);   
   IntegerVector O=bincounter_cpp(x, bins);
   NumericVector p(bins.size());
   Rcpp::Environment base("package:base");
   Rcpp::Function formals_r = base["formals"];
   Rcpp::List res = formals_r(Rcpp::_["fun"]=pnull);
   if(res.size()==1) p=pnull(bins);
   else p=pnull(bins, param);
   if(rate==0) rate=n;
   for(i=0;i<k;++i) E[i] = rate*(p[i+1]-p[i])/(p[k]-p[0]);
   double chi = 0;
   for(i=0;i<E.size();++i) {
     if(formula=="Pearson") chi = chi + (O[i]-E[i])*(O[i]-E[i])/E[i];
     else 
       if(O[i]>0) chi = chi + 2.0*(E[i]-O[i]+O[i]*log(O[i]/E[i]));
          
   }  
   return chi;
}

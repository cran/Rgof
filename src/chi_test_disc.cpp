#include <Rcpp.h>
#include "binner_disc.h"
#include "chi_stat_disc.h"

using namespace Rcpp;
//'
//' This function performs a number of chisquare gof tests for discrete data
//' @param  x data set
//' @param  pnull  cdf under the null hypothesis
//' @param  param  starting values of multi-D minimum chi square minimization
//' @param  nbins =c(100, 10) number of bins for chisquare tests
//' @param  formula type of chi square formula to use
//' @param  rate rate of Poisson if sample size is random, 0 otherwise
//' @param  Minimize Should minimum chi square be found?
//' @param minexpcount =2 minimal expected bin count required
//' @keywords internal
//' @return A numeric matrix of test statistics, degrees of freedom and p values
// [[Rcpp::export]]
Rcpp::NumericMatrix chi_test_disc (
    Rcpp::IntegerVector x, 
    Rcpp::Function pnull, 
    Rcpp::NumericVector param,
    Rcpp::IntegerVector nbins= Rcpp::IntegerVector::create(100, 10),
    std::string formula="Pearson",
    double rate=0.0,
    int Minimize=0,
    double minexpcount=2.0) {
  int i, nb;
  Rcpp::NumericMatrix out(2, 3);
  Rcpp::NumericVector p(x.size());
  Rcpp::List res;
  colnames(out) = Rcpp::CharacterVector::create("Statistic", "df", "p-value");
  rownames(out) = Rcpp::CharacterVector::create("chi large", "chi small");
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  res = formals_r(Rcpp::_["fun"]=pnull);
  if(res.size()==0) p=pnull();
  else p=pnull(param);
  for(i=0;i<=1;++i) {
    nb = nbins[i];
    if(nb>x.size()) nb=x.size();
    IntegerVector bins = binner_disc(x, p, nb, minexpcount);
    if(Minimize==0) {
       out(i, 0) = chi_stat_disc(param, x, pnull, bins, formula, rate);         
       out(i, 1) = bins.size()-1;
       out(i, 2) = 1-R::pchisq(out(i, 0), out(i, 1), 1, 0);
    } else {
        Rcpp::Environment base("package:stats");
        Rcpp::Function optim_r = base["optim"]; 
        res = optim_r(
            Rcpp::_["par"]=param, 
            Rcpp::_["fn"]=Rcpp::InternalFunction(&chi_stat_disc), 
            Rcpp::_["method"]="BFGS",
            Rcpp::_["x"]=x, 
            Rcpp::_["pnull"]=pnull, 
            Rcpp::_["bins"]=bins,
            Rcpp::_["formula"]=formula,
            Rcpp::_["rate"]=rate);
        out(i, 0) = res[1];
        if(rate==0) out(i, 1) = bins.size()-1-param.size();
        else out(i, 1) = bins.size()-param.size();
        out(i, 2) = 1-R::pchisq(out(i, 0), out(i, 1), 1, 0);
    }
   }   
   return out;
  
}

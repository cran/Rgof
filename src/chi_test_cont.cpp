#include <Rcpp.h>
#include "binner_cont.h"
#include "chi_stat_cont.h"
#include <string>
using namespace std;
using namespace Rcpp;
//'
//' This function performs a number of chisquare gof tests for continuous data
//' @param  x data set
//' @param  pnull  cdf under the null hypothesis
//' @param  param  starting values of multi-D minimum chi square minimization
//' @param  nbins =c(100, 10) number of bins for chisquare tests
//' @param  formula Formula of chi square to use
//' @param  rate rate of Poisson if sample size is random
//' @param  Range  =(-99999, 99999) limits of possible observations, if any
//' @param  Minimize Should minimum chi square be found?
//' @keywords internal
//' @return A numeric matrix of test statistics, degrees of freedom and p values
// [[Rcpp::export]]
Rcpp::NumericMatrix chi_test_cont (
    Rcpp::NumericVector x, 
    Rcpp::Function pnull, 
    Rcpp::NumericVector param,
    std::string formula="Pearson",
    double rate=0.0,
    Rcpp::IntegerVector nbins= Rcpp::IntegerVector::create(100, 10),
    Rcpp::NumericVector Range=Rcpp::NumericVector::create(-99999, 99999),
    int Minimize=0
  ) {
  int i, j, nb, which;
  Rcpp::List res;
  Rcpp::NumericMatrix out(4, 3);
  colnames(out) = Rcpp::CharacterVector::create("Statistic", "df", "p-value");
  rownames(out) = Rcpp::CharacterVector::create("EP large", "ES large", "EP small", "ES small");
 
  std::sort(x.begin(), x.end());  
  j = 0;
  for(i=0;i<=1;++i) {
     nb = nbins[i];
     if(nb>x.size()/5.0) nb=x.size()/5.0;
     for(which=1;which<=2;++which) {
       NumericVector bins = binner_cont(x, pnull, param, nb, which, Range);

       if(Minimize==0) {
          out(j, 0) = chi_stat_cont(param, x, pnull, bins, formula, rate);
          out(j, 1) = bins.size()-2;
          out(j, 2) = 1-R::pchisq(out(j, 0), out(j, 1), 1, 0);
       } 
       else {
         Rcpp::Environment base("package:stats");
         Rcpp::Function optim_r = base["optim"];           
         res = optim_r(
                Rcpp::_["par"]=param, 
                Rcpp::_["fn"]=Rcpp::InternalFunction(&chi_stat_cont), 
                Rcpp::_["method"]="BFGS",
                Rcpp::_["x"]=x, 
                Rcpp::_["pnull"]=pnull, 
                Rcpp::_["bins"]=bins,
                Rcpp::_["formula"]=formula, 
                Rcpp::_["rate"]=rate);
         out(j, 0) = res[1];
         if(rate==0) out(j, 1) = bins.size()-2-param.size();
         else out(j, 1) = bins.size()-1-param.size();
         out(j, 2) = 1-R::pchisq(out(j, 0), out(j, 1), 1, 0);
       }
       j = j+1;
     }
     
   }   
   return out;
  
}

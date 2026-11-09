#include <Rcpp.h>
#include "calcTS.h"
using namespace Rcpp;

//' helper functions to do p value adjustment
//' @param dta data set
//' @param rnull R function (generate data under null hypothesis)
//' @param vals values of discrete random variable
//' @param TS function to calculate test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param B  =1000 Number of simulation runs
//' @keywords internal
//' @return A matrix of powers
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix gof_adj_C1(
        List dta,
        Rcpp::Function rnull,
        Rcpp::NumericVector vals,
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,
        const int B=1000) {

  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
  NumericVector x=dta["x"];
  Function phat=TSextra["phat"];
  NumericVector p=phat(x);
  NumericVector TS_sim=calcTS(dta, TS, typeTS, TSextra);
  int  i, num_tests=TS_sim.size();
  NumericVector xsim(num_tests);
  List dtasim=List::create(Named("x")=xsim, Named("vals")=vals);
  NumericMatrix A(B, num_tests);
  for(i=0;i<B;++i) {
      if(resrnull.size()==0) xsim=rnull();
      else xsim=rnull(p);
      dtasim["x"]=xsim;
      TS_sim=calcTS(dtasim, TS, typeTS, TSextra);
      A(i, _)=TS_sim;
  }
  
  return A;
}

//' helper functions to do p value adjustment
//' @param dta data set
//' @param rnull R function (generate data under null hypothesis)
//' @param vals values of discrete random variable
//' @param TS function to calculate test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param A matrix of test statistic values 
//' @param B  =1000 Number of simulation runs
//' @keywords internal
//' @return A matrix of powers
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix gof_adj_C2(
     List dta,
     Rcpp::Function rnull,
     Rcpp::NumericVector vals,
     Rcpp::Function TS,
     int typeTS,
     Rcpp::List TSextra,
     Rcpp::NumericMatrix A,
     const int B=1000) {
   
   Rcpp::Environment base("package:base");
   Rcpp::Function formals_r = base["formals"];
   Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
   NumericVector x=dta["x"];
   Function phat=TSextra["phat"];
   NumericVector p=phat(x);
   NumericVector TS_sim=calcTS(dta, TS, typeTS, TSextra);
   int  i, j, k, m=A.nrow(), num_tests=TS_sim.size();
   NumericVector xsim(num_tests);
   double dd=1.0/double(m);
   List dtasim=List::create(Named("x")=xsim, Named("vals")=vals);
   NumericMatrix pvals(B+1, num_tests);
/* row 0: the observed data, rows 1,..,B: data simulated under H0. 
   p value = proportion of the values in A at least as large as the 
   test statistic. */
   for(i=0;i<=B;++i) {
       if(i==0) TS_sim=calcTS(dta, TS, typeTS, TSextra);
       else {
         if(resrnull.size()==0) xsim=rnull();
         else xsim=rnull(p);
         dtasim["x"]=xsim;
         TS_sim=calcTS(dtasim, TS, typeTS, TSextra);
       }      
       for(j=0;j<num_tests;++j) {       
          pvals(i,j)=0.0;   
          double tol=1e-9*std::abs(TS_sim(j));
          for(k=0;k<m;++k) {
             if(A(k,j)>=TS_sim(j)-tol)  pvals(i, j)=pvals(i, j)+dd;
         }
       }   
   }
   return pvals;
 }

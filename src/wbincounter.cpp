#include <Rcpp.h>
#include "Cpporder.h"
using namespace Rcpp;

//' Find counts or sum of weights in bins. Useful for power calculations. Replaces hist command from R.
//' 
//' @param x numeric vector
//' @param bins numeric vector
//' @param w numeric vector of weights 
//' @keywords internal
//' @return sum of weights in bins
// [[Rcpp::export]]
Rcpp::NumericVector wbincounter(Rcpp::NumericVector x, 
                                Rcpp::NumericVector bins, 
                                Rcpp::NumericVector w) {
  int n=x.size(), m=bins.size(), i, j;
  NumericVector wc(m-1);
  
  w = Cpporder(w, x);
  x = Cpporder(x, x);
  i=0;
  j=0;
  while ( (j<m-1) && (i<n) ) {
    if(x(i)<=bins(j+1)) {
       wc(j) = wc(j) + w(i);
       ++i;
    }  
    else ++j;
    
  }
  return wc;   
}

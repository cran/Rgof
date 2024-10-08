#include <Rcpp.h>
#include "Cpporder.h"
using namespace Rcpp;

//' count events in bins. Useful for power calculations. Replaces hist command from R.
//' 
//' @param x numeric vector
//' @param bins numeric vector
//' @keywords internal
//' @return Integer vector of counts
// [[Rcpp::export]]
Rcpp::IntegerVector bincounter(Rcpp::NumericVector x, Rcpp::NumericVector bins) {
  int n=x.size(), m=bins.size(), i, j;
  Rcpp::IntegerVector xc(m-1);
  x = Cpporder(x, x);

  i=0;
  j=0;
  while ( (j<m-1) && (i<n) ) {
    if(x(i)<=bins(j+1)) {
      ++xc(j);
      ++i;
    }  
    else ++j;
    
  }

  return xc;   
}

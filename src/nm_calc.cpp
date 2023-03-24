#include <Rcpp.h>

using namespace Rcpp;

//' calculate a matrix of numbers needed for Zhangs tests for discrete data
//' 
//' @param n An integer.
//' @keywords internal
//' @return A matrix of numbers
// [[Rcpp::export]]
NumericMatrix nm_calc(int n) {
  double tmp;
  NumericMatrix nm(n+1, 4);

  nm(0, 0) = 0.0;
  nm(0, 1) = 0.0;
  nm(0, 2) = 0.0;
  nm(0, 3) = 0.0;
  
  nm(1, 0) = 1.0/(n-0.5);
  nm(1, 1) = 1.0/0.5;
  tmp = log( (n-0.5)/(0.25)-1 );
  nm(1, 2) = tmp;
  nm(1, 3) = tmp*tmp;
  for(int j=1;j<n;++j) {
    nm(j+1, 0) = nm(j, 0) + 1.0/(n-double(j)-0.5);
    nm(j+1, 1) = nm(j, 1) + 1.0/(double(j)+0.5);
    tmp = log( (n-0.5)/(double(j)+0.25)-1 );
    nm(j+1, 2) = nm(j, 2) + tmp;
    nm(j+1, 3) = nm(j, 3) + tmp*tmp;
  }

  return nm;
}

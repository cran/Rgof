#include <Rcpp.h>
using namespace Rcpp;

//' sort vector y by values in vector x
//' 
//' @param y numeric vector
//' @param x numeric vector
//' @keywords internal
//' @return numeric vector
// [[Rcpp::export]]
NumericVector Cpporder(NumericVector y, NumericVector x) {
  std::vector<int> Index(x.size());
  std::iota(Index.begin(), Index.end(), 0);
  std::sort(Index.begin(), Index.end(),
            [&](int A, int B) -> bool {
              return x[A] < x[B];
            });
  NumericVector out(x.size());
  for(int i=0;i<x.size();++i) out[i]=y[Index[i]];
  return out;
}

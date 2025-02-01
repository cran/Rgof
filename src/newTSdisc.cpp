#include <Rcpp.h>
using namespace Rcpp;

//' a local function needed for the vignette
//' 
//' @param x An integer vector.
//' @param pnull cdf.
//' @param param parameters for pnull in case of parameter estimation.
//' @param vals A numeric vector with the values of the discrete rv.
//' @return A vector with test statistics
//' @export
// [[Rcpp::export]]
NumericVector newTSdisc(IntegerVector x, 
                      Function pnull,  
                      NumericVector param,
                      NumericVector vals) {
    
  Rcpp::CharacterVector methods=CharacterVector::create("CvM alt");    
  int const nummethods=methods.size();
  int k=x.size(), n, i;
  NumericVector TS(nummethods), ecdf(k), Fx(k);
  double tmp;
  TS.names() =  methods;
  Fx=pnull(param);
  n=0;
  for(i=0;i<k;++i) n = n + x[i];
  ecdf(0) = double(x(0))/double(n);  
  for(i=1;i<k;++i) {
    ecdf(i) = ecdf(i-1) + x(i)/double(n);
  }

  tmp = std::abs(ecdf[0]-Fx(0))*Fx(0);
  for(i=1;i<k;++i) 
     tmp = tmp + std::abs(ecdf(i)-Fx(i))*(Fx(i)-Fx(i-1));
  TS(0) = tmp;
 
  return TS;
}

#include <Rcpp.h>
using namespace Rcpp;

//' routine to bin discrete data and assure expected counts>5 
//' 
//' @param  x An integer vector of counts.
//' @param  p cumulative distribution function values at vals
//' @param  k number of desired bins
//' @param  minexpcount =2 minimal expected bin count required
//' @keywords internal
//' @return An integer vector of indices
// [[Rcpp::export]]
Rcpp::IntegerVector binner_disc(
    Rcpp::IntegerVector x, 
    Rcpp::NumericVector p, 
    int k=10,
    double minexpcount=2) {
  int n=sum(x), i;
  double tmp;
  NumericVector E(x.size());  
  IntegerVector out(x.size());
  for(i=0;i<x.size();++i) out[i]=i;
  E[0] = double(n)*p[0];
  for(i=1;i<x.size();++i) E[i] = double(n)*(p[i]-p[i-1]);
  E = E/max(p);
  
  if(k>x.size()) k = x.size();
  
  while ( (E.size()>k) || (min(E)<minexpcount)) {
      int whichmin=0;
      tmp=E[0];
      for(i=1;i<E.size();++i) {
         if(tmp>E[i]) {
           tmp=E[i];
           whichmin=i;
         }
      }
   
      if(whichmin<E.size()/2.0) {
        E[whichmin]=E[whichmin]+E[whichmin+1];
        E.erase(whichmin+1);
        out.erase(whichmin+1);
      }
      else {
        E[whichmin]=E[whichmin]+E[whichmin-1];
        E.erase(whichmin-1);
        out.erase(whichmin-1);
      }

  }
  return out;

}

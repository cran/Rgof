#include <Rcpp.h>
using namespace Rcpp;

//' routine to bin continuous data and assure expected counts>5 
//' 
//' @param  x A numeric vector of numbers.
//' @param  pnull function to find CDF
//' @param  param parameter of pnull
//' @param  k =10, number of bins
//' @param  which =1, 1=equal probability bins, 2-equal size bins
//' @param Range =(-99999,99999) limits of possible observations, if any
//' @keywords internal
//' @return A vectors of bin edges
// [[Rcpp::export]]
Rcpp::NumericVector binner_cont(
    Rcpp::NumericVector x, 
    Rcpp::Function pnull, 
    Rcpp::NumericVector param, 
    int k=10,
    int which=1,
    Rcpp::NumericVector Range=Rcpp::NumericVector::create(-99999,99999)) {
  const int minE=2;  
  int n=x.size(), i;
  NumericVector bins(k+1), E(k);  
  if(Range[0]!=-99999) bins[0]=Range[0];
  else bins[0]=x[0]-1e-10;
  if(Range[1]!=99999) bins[k]=Range[1];
  else bins[k]=x[x.size()-1]+1e-10;
  if(which==1)
    for(i=1;i<k;++i) bins[i] = x[(n-1)*i/k]+1e-10;
  else
    for(i=1;i<k;++i) bins[i] = bins[0]+double(i)/k*(bins[k]-bins[0]);

  
  NumericVector p(bins.size());
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List res = formals_r(Rcpp::_["fun"]=pnull);
  if(res.size()==1) p=pnull(bins);
  else p=pnull(bins, param);
  for(i=0;i<k;++i) E[i] = n*(p[i+1]-p[i])/(p[k]-p[0]);   
  
  while ( min(E)<minE ) {
    int whichmin=0;
    double tmp=E[0];
    for(i=1;i<E.size();++i) {
      if(tmp>E[i]) {
        tmp=E[i];
        whichmin=i;
      }
    }
    if(whichmin>0) {
      E[whichmin]=E[whichmin-1]+E[whichmin];
      E.erase(whichmin-1);
      bins.erase(whichmin);
    }
    else {
      E[0]=E[0]+E[1];
      E.erase(1);
      bins.erase(1);
    }
  }  
  return bins;
}

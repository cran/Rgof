#include <Rcpp.h>
using namespace Rcpp;
//' This function calculates the test statistics
//' @param  dta list with data set
//' @param  pnull  cdf under the null hypothesis
//' @param  param estimated parameters (or some constant)
//' @param  TS routine
//' @param  typeTS format of TS 
//' @param  TSextra list passed to TS function
//' @keywords internal
//' @return A vector of numbers
// [[Rcpp::export]]
NumericVector calcTS(
     Rcpp::List dta,
     Rcpp::Function pnull,
     Rcpp::NumericVector param,
     Rcpp::Function TS,      
     int typeTS,
     Rcpp::List TSextra) {
   NumericVector TS_data, wx;
   if(typeTS==1) {
     Function qnull=TSextra["qnull"]; 
     TS_data=TS(dta["x"], pnull, param, qnull);
   }   
   if(typeTS==2) {
     Function w=TSextra["w"];
     if(std::abs(param(0)+99)<0.01) wx=w(dta["x"]);
     else wx=w(dta["x"], param);
     TS_data=TS(dta["x"], pnull, param, wx);
   }  
   if(typeTS==3) TS_data=TS(dta["x"], pnull, param);
   if(typeTS==4) TS_data=TS(dta["x"], pnull, param, TSextra);
   if(typeTS==5) TS_data=TS(dta["x"], pnull, param, dta["vals"]); 
   if(typeTS==6) TS_data=TS(dta["x"], pnull, param, dta["vals"], TSextra);
   return  TS_data;
 }
 


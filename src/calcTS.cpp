#include <Rcpp.h>
using namespace Rcpp;
//' This function calculates the test statistics for  data
//' @param  dta data set as a list
//' @param  TS routine
//' @param  typeTS format of TS
//' @param  TSextra list passed to TS function
//' @keywords internal
//' @return A vector of numbers
// [[Rcpp::export]]
NumericVector calcTS(
     Rcpp::List dta, 
     Rcpp::Function TS,
     int typeTS,
     Rcpp::List TSextra) {
   NumericVector TS_data;
   Function phat=TSextra["phat"];
   Function w=TSextra["w"];
   Rcpp::Environment base("package:base");
   Rcpp::Function formals_r = base["formals"];
   Rcpp::List resw = formals_r(Rcpp::_["fun"]=w);
   NumericVector p=phat(dta["x"]);
   NumericVector wx;
   if(typeTS==2) {
     if(resw.size()==1) wx=w(dta["x"]);
     else wx=w(dta["x"],p);
   }   
   if(typeTS<5) {
     if(typeTS==1) TS_data=TS(dta["x"], TSextra["pnull"], p, TSextra["qnull"]); 
     if(typeTS==2) TS_data=TS(dta["x"], TSextra["pnull"], p, wx); 
     if(typeTS==3) TS_data=TS(dta["x"], TSextra["pnull"], p);   
     if(typeTS==4) TS_data=TS(dta["x"], TSextra["pnull"], p, TSextra);
   }
   else {
     if(typeTS==5) TS_data=TS(dta["x"], TSextra["pnull"], p, dta["vals"]);
     if(typeTS==6) TS_data=TS(dta["x"], TSextra["pnull"], p, dta["vals"], TSextra);  
   }
   return  TS_data;
 }
 
 

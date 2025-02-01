#ifndef NEWTSDISC_H
#define NEWTSDISC_H

#include <Rcpp.h>
Rcpp::NumericVector newTSdisc(Rcpp::IntegerVector x, 
                       Rcpp::Function pnull, 
                       Rcpp::NumericVector param,
                       Rcpp::NumericVector vals                        
                   );

#endif

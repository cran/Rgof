#' This function finds the power of various chi-square tests for continuous data
#' @param  pnull function to find cdf under  null hypothesis
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  qnull =NA function to find quantiles under null hypothesis, if available
#' @param  phat =function(x) -99, function to estimate parameters 
#' @param  w =function(x) -99, optional weight function
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Range  =c(-99999, 99999) limits of possible observations, if any
#' @param  B =1000 number of simulation runs to find power
#' @param  nbins =c(50,10), number of bins for chi square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat =TRUE, if TRUE param is estimated parameters and no minimization is used
#' @return A numeric matrix of power values.

chi_power_cont = function(
        pnull, ralt, param_alt, qnull=NA,
        phat=function(x) -99, w=function(x) -99,
        alpha=0.05, Range=c(-99999, 99999), B=1000,  
        nbins=c(50, 10), rate=0, minexpcount=5, ChiUsePhat=TRUE) {

   pwr=matrix(0, length(param_alt), 8)
   colnames(pwr) = c("ES-l-P", "ES-s-P", "EP-l-P", "EP-s-P",
                     "ES-l-L", "ES-s-L", "EP-l-L", "EP-s-L")
   rownames(pwr) = param_alt
   for(i in 1:length(param_alt)) {
     for(j in 1:B) {
         x = ralt(param_alt[i])
         tmp = c(chi_test_cont(x, pnull, w, phat, qnull, 
           nbins, rate, Range, minexpcount, ChiUsePhat)[, 2])
         pwr[i, ] = pwr[i, ] + ifelse(tmp<alpha, 1, 0)
     }
     pwr[i, ] = pwr[i, ]/B
   }
   pwr
}

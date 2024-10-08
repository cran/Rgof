#' This function finds the power of various chi-square tests for continuous data
#' @param  pnull function to find cdf under  null hypothesis
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat =function(x) -99, routine to estimate parameters
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B =1000 number of simulation runs to find power
#' @param  nbins =c(50,10), number of bins for chi square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, should chi square use minimum chi square method?
#' @return A numeric matrix of power values.

chi_power_disc = function(pnull, ralt, param_alt, 
                phat=function(x) -99, alpha=0.05, B=1000,  
                nbins=c(50, 10), rate=0, minexpcount=5, 
                ChiUsePhat=TRUE) {
   pwr=matrix(0, length(param_alt), 4)
   colnames(pwr) = c(c("l-P", "s-P", "l-L", "s-L"))
   rownames(pwr) = param_alt
   for(i in 1:length(param_alt)) {
     x = ralt(param_alt[i])
     allbins = make_bins_disc(x, pnull, phat, nbins=nbins, minexpcount=minexpcount)
     for(j in 1:B) {
         x = ralt(param_alt[i])
         tmp = chi_test_disc(x, pnull, phat, 
                  ChiUsePhat=ChiUsePhat, allbins=allbins)[, 2]
         pwr[i, ] = pwr[i, ] + ifelse(tmp<alpha, 1, 0)
     }
     pwr[i, ] = pwr[i, ]/B
   }
   pwr
}

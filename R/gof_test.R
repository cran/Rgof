#' This function performs a number of gof tests
#' @param  x data set
#' @param  vals =NA, values of discrete RV, or NA if data is continuous
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  w (Optional) function to calculate weights, returns -99 if no weights
#' @param  phat  =function(x) -99, function to estimate parameters from the data, or -99 if no parameters are estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  TSextra =NA, list passed to TS, if desired, or NA
#' @param  nbins =c(100, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any, for chi-square tests
#' @param  B   =5000  number of simulation runs
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  maxProcessors =1, number of processors to use in parallel processing. 
#' @param  doMethods Methods to include in tests
#' @return A list with vectors of test statistics and p.values
#' @export
#' @examples
#' # Tests to see whether data comes from a standard normal distribution.
#' pnull = function(x) pnorm(x)
#' rnull = function()  rnorm(100)
#' x = rnorm(100)
#' gof_test(x, NA, pnull, rnull)
#' # Tests to see whether data comes from a normal distribution with standard deviation 1 
#' # and the mean estimated.
#' pnull=function(x, m) pnorm(x, m)
#' rnull=function(m) rnorm(100, m)
#' TSextra = list(qnull=function(x, m=0) qnorm(x, m), 
#'           pnull=function(x, m=0) pnorm(x, m), phat=function(x) mean(x))
#' phat=function(x) mean(x)
#' x = rnorm(100, 1, 2)
#' gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra)
#' # Tests to see whether data comes from a binomial (10, 0.5) distribution.
#' vals=0:10
#' pnull = function() pbinom(0:10, 10, 0.5)
#' rnull = function() table(c(0:10, rbinom(1000, 10, 0.5)))-1
#' x = rnull() 
#' gof_test(x, vals, pnull, rnull, doMethods="all")
#' # Tests to see whether data comes from a binomial distribution with 
#' # the success probability estimated from the data.
#' pnull = function(p=0.5) pbinom(0:10, 10, ifelse(p>0&&p<1, p, 0.001))
#' rnull = function(p=0.5) table(c(0:10, rbinom(1000, 10, 
#'                   ifelse(p>0&&p<1, p, 0.001))))-1
#' phat=function(x) mean(rep(0:10,x))/10 
#' gof_test(x, vals, pnull, rnull, phat=phat) 
#'
gof_test <- function(x, vals= NA, pnull, rnull, 
                    w=function(x) -99, phat=function(x) -99, 
                    TS, TSextra=NA, nbins=c(50, 10), rate=0, 
                    Range=c(-Inf, Inf), B=5000,  minexpcount=5.0,  
                    ChiUsePhat=TRUE, maxProcessors=1, doMethods="all") {

  if(length(x)>10000 && maxProcessors==1)
      message("Consider using parallel processing with maxProcessor= (your number of cores)")
# adjust number of bins to account for parameter estimation
   if(abs(phat(x)[1]+99)<0.001) nbins=nbins+length(phat(x)) 
   if(any(is.na(vals))) { # continuous data/model 
     # do some checks to see arguments are given correctly
     check.functions(pnull, rnull, phat, x=x)
     if(missing(TS)) # run built-in methods
        out = gof_test_cont(x, pnull, rnull, w, phat, TSextra=TSextra, nbins=nbins, 
                  rate=rate, Range=Range, B=B, minexpcount=minexpcount, ChiUsePhat=ChiUsePhat, 
                  maxProcessors=maxProcessors, doMethods=doMethods)
     else # run user-provided tests
       out = gof_test_cont(x, pnull, rnull, w, phat, TS=TS, TSextra=TSextra, nbins=nbins, 
                    rate=rate, Range=Range, B=B, minexpcount=minexpcount, 
                    ChiUsePhat=ChiUsePhat, maxProcessors=maxProcessors, doMethods=doMethods)           
   }
   else { #discrete data/model
     # do some checks to see arguments are given correctly
     check.functions(pnull, rnull, phat, vals, x)
     if(missing(TS)) # run built-in methods
     out = gof_test_disc(x, pnull, rnull, vals, phat, TSextra=TSextra, nbins=nbins, 
                   rate=rate, B=B, minexpcount=minexpcount, ChiUsePhat=ChiUsePhat,
                   maxProcessors=maxProcessors, doMethods=doMethods)     
     else # run user-provided tests
       out = gof_test_disc(x, pnull, rnull, vals, phat, TS=TS, TSextra=TSextra, nbins=nbins, 
                   rate=rate, B=B, minexpcount=minexpcount, ChiUsePhat=ChiUsePhat,
                   maxProcessors=maxProcessors, doMethods=doMethods)     
   }
   # make output look nice
   signif.digits(out)
}

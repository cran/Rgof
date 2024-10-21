#' Find the power of various gof tests for continuous data.
#' @param  pnull function to find cdf under  null hypothesis
#' @param  vals =NA, values of rv, if data is discrete, NA if data is continuous
#' @param  rnull function to generate data under  null hypothesis
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  w (Optional) function to calculate weights, returns -99 if no weights
#' @param  phat =function(x) -99 function to estimate parameters from the data, or -99
#' @param  TS user supplied function to find test statistics
#' @param  TSextra =NA, list provided to TS
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any
#' @param  B =c(1000, 1000), number of simulation runs to find power and null distribution
#' @param  nbins =c(100,10), number of bins for chi square tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  maxProcessors maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @return A numeric matrix of power values.
#' @export 
#' @examples
#' # Power of tests when null hypothesis specifies the standard normal distribution but 
#' # true data comes from a normal distribution with mean different from 0.
#' pnull = function(x) pnorm(x)
#' rnull = function()  rnorm(50)
#' ralt = function(mu)  rnorm(50, mu)
#' TSextra = list(qnull=function(x) qnorm(x))
#' gof_power(pnull, NA, rnull, ralt, c(0.25, 0.5), TSextra=TSextra, B=c(500, 500))
#' # Power of tests when null hypothesis specifies normal distribution and 
#' # mean and standard deviation are estimated from the data. 
#' # Example is not run because it takes several minutes.
#' # true data comes from a normal distribution with mean different from 0.
#' pnull = function(x, p=c(0, 1)) pnorm(x, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' rnull = function(p=c(0, 1))  rnorm(50, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' phat = function(x) c(mean(x), sd(x))
#' TSextra = list(qnull = function(x, p=c(0, 1)) qnorm(x, p[1],  
#'                ifelse(p[2]>0.001, p[2], 0.001))) 
#' \donttest{gof_power(pnull, NA, rnull, ralt, c(0, 1), phat=phat, TSextra=TSextra, 
#'           B=c(200, 200), maxProcessor=2)}
#' # Power of tests when null hypothesis specifies Poisson rv with rate 100 and 
#' # true rate is 100.5
#' vals = 0:250
#' pnull = function() ppois(0:250, 100)
#' rnull =function () table(c(0:250, rpois(1000, 100)))-1
#' ralt =function (p) table(c(0:250, rpois(1000, p)))-1
#' gof_power(pnull, vals, rnull, ralt, param_alt=100.5,  B=c(500,500))
#' # Power of tests when null hypothesis specifies a Binomial n=10 distribution 
#' # with the success probability estimated
#' vals = 0:10
#' pnull=function(p) pbinom(0:10, 10, ifelse(0<p&p<1, p, 0.001))
#' rnull=function(p) table(c(0:10, rbinom(1000, 10, ifelse(0<p&p<1, p, 0.001))))-1
#' ralt=function(p) table(c(0:10, rbinom(1000, 10, p)))-1
#' phat=function(x) mean(rep(0:10,x))/10
#' \donttest{gof_power(pnull, vals, rnull, ralt, c(0.5, 0.6), phat=phat,
#'                     B=c(200, 200), maxProcessor=2)}
#'
gof_power=function(pnull, vals=NA, rnull, ralt, param_alt, 
        w=function(x) -99, phat=function(x) -99, TS, TSextra=NA, 
        alpha=0.05, Range  =c(-Inf, Inf), B=c(1000, 1000),nbins=c(50,10), 
        rate=0, maxProcessors, minexpcount=5.0, ChiUsePhat=TRUE) {
  
  x = ralt(param_alt[1]) # get an example data set
  # adjust number of bins to account for parameter estimation
  if(abs(phat(x)[1]+99)<0.001) nbins=nbins+length(phat(x)) 
  if(any(is.na(vals))) { # continuous data/model
    check.functions(pnull, rnull, phat, x=x) # do some sanity checks
    if(missing(TS)) # use built-in tests
        out = gof_power_cont(pnull, rnull, ralt, param_alt, w, phat,  
                         TSextra=TSextra, alpha=alpha, Range=Range, B=B, 
                         nbins=nbins, rate=rate, maxProcessors=maxProcessors, 
                         minexpcount=minexpcount, ChiUsePhat=ChiUsePhat)
    else # do user-provided tests
       out = gof_power_cont(pnull, rnull, ralt, param_alt, w, phat,  TS=TS,
                           TSextra=TSextra, alpha=alpha, Range=Range, B=B, 
                           nbins=nbins, rate=rate, maxProcessors=maxProcessors, 
                           minexpcount=minexpcount, ChiUsePhat=ChiUsePhat)  
  }
  else { # discrete data
    check.functions(pnull, rnull, vals=vals, phat=phat, x=x)
    if(missing(TS))
      out = gof_power_disc(pnull, rnull, vals, ralt, param_alt, phat,  
                           TSextra=TSextra, alpha=alpha, B=B, 
                           nbins=nbins, rate=rate, 
                           maxProcessors=maxProcessors, 
                           minexpcount=minexpcount,
                           ChiUsePhat=ChiUsePhat)
    else
      out = gof_power_disc(pnull, rnull, vals, ralt, param_alt, phat,  TS=TS,
                           TSextra=TSextra, alpha=alpha, B=B, 
                           nbins=nbins, rate=rate, 
                           maxProcessors=maxProcessors, 
                           minexpcount=minexpcount,
                           ChiUsePhat=ChiUsePhat)     
  }
  out
}

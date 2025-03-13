#' Find the power of various gof tests for discrete data.
#' @param  pnull cumulative distribution function under the null hypothesis
#' @param  rnull  a function to generate data under  null hypothesis
#' @param  vals values of discrete rv.
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat  =function(x) -99, function to estimate parameters from the data, -99 if no parameters are estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  TSextra =NA, list passed to TS, if desired
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B =1000 number of simulation runs
#' @param  nbins =c(100, 10) number of bins for chi square tests
#' @param  rate  rate of Poisson if sample size is random
#' @param  maxProcessor maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  minexpcount =5 minimal number of expected counts in each bin for chi square tests
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @return A numeric matrix of power values.

gof_power_disc=function(pnull, rnull, vals, ralt, param_alt, phat=function(x) -99, 
        TS, TSextra=NA, alpha=0.05, B=1000,  nbins=c(100,10), 
        rate=0, maxProcessor, minexpcount=5.0, ChiUsePhat=TRUE) {

    x = ralt(param_alt[1])
    WithEstimation = TRUE
    if(phat(x)[1]==-99) WithEstimation = FALSE

    if(any(is.na(TSextra))) TSextra = list(p=phat(x))
    else TSextra = c(TSextra, p=phat)
    
    if(missing(TS)) {
      typeTS=5
      TS =  TS_disc
      TS_data = TS(x, pnull, phat(x), vals)
    }  
    else {
      if(length(formals(TS))==4) {
        typeTS=5
        TS_data = TS(x, pnull, phat(x), vals)
      }
      if(length(formals(TS))==5) {
        typeTS=6
        TS_data = TS(x, pnull, phat(x), vals, TSextra)
      }
      if(length(formals(TS))>6) {
        message("TS should have either 4 or 5 arguments")
        return(NULL)
      } 
      if(is.null(names(TS_data))) {
        message("result of TS has to be a named vector")
        return(NULL)
      }
      if(substr(deparse(TS)[2], 1, 5)==".Call") {
        message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
        maxProcessor=1
      }  
    }
    nummethods = length(TS_data)
    methods = names(TS_data)
    # Do chi square tests if built-in TS is used  
    chiout=NULL
    if(typeTS==0) { #Run chi square tests
      chiout = chi_power_disc(pnull, ralt, param_alt, 
                              phat, alpha , B, 
                              nbins, rate, minexpcount,
                              ChiUsePhat)[,1:2, drop=FALSE]
    }
# Now the other tests    
    out = power_disc_R(pnull = pnull,
               rnull = rnull, 
               vals = vals,
               ralt = ralt, 
               param_alt = param_alt, 
               phat = phat,
               TS = TS,
               typeTS = typeTS,
               TSextra = TSextra,
               B = B,
               alpha = alpha,
               maxProcessor = maxProcessor)
    rownames(out) = param_alt 
    colnames(out) = methods
    out = cbind(out, chiout)
    if(is.matrix(out) && length(param_alt)==1) out=out[1, ]
    out
}

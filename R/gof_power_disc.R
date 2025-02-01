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
#' @param  B =c(1000, 1000), number of simulation runs to find power and null distribution
#' @param  nbins =c(100, 10) number of bins for chi square tests
#' @param  rate  rate of Poisson if sample size is random
#' @param  maxProcessors maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  minexpcount =5 minimal number of expected counts in each bin for chi square tests
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @return A numeric matrix of power values.

gof_power_disc=function(pnull, rnull, vals, ralt, param_alt, phat=function(x) -99, 
        TS, TSextra=NA,
        alpha=0.05, B=c(1000, 1000),  nbins=c(100,10), 
        rate=0, maxProcessors, minexpcount=5.0, ChiUsePhat=TRUE) {

    x = ralt(param_alt[1])
    WithEstimation = TRUE
    if(phat(x)[1]==-99) WithEstimation = FALSE

    if(any(is.na(TSextra))) TSextra = list(p=phat(x))
    else TSextra = c(TSextra, p=phat)
    
    if(missing(TS)) {
      typeTS=0
      TS =  TS_disc
      TS_data = TS(x, pnull, phat(x), vals)
    }  
    else {
      if(length(formals(TS))==4) {
        typeTS=1
        TS_data = TS(x, pnull, phat(x), vals)
      }
      if(length(formals(TS))==5) {
        typeTS=2
        TS_data = TS(x, pnull, phat(x), vals, TSextra)
      }
      if(length(formals(TS))>5) {
        message("TS should have either 4 or 5 arguments")
        return(NULL)
      } 
      if(is.null(names(TS_data))) {
        message("result of TS has to be a named vector")
        return(NULL)
      }
      if(substr(deparse(TS)[2], 1, 5)==".Call") {
        message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
        maxProcessors=1
      }  
    }
    nummethods = length(TS_data)
    methods = names(TS_data)
    # Do chi square tests if built-in TS is used  
    chiout=NULL
    if(typeTS==0) { #Run chi square tests
      chiout = chi_power_disc(pnull, ralt, param_alt, 
                              phat, alpha , B[2], 
                              nbins, rate, minexpcount,
                              ChiUsePhat)[,1:2, drop=FALSE]
    }
# Now the other tests    
# With parameter estimation    
    if(WithEstimation) {
      if(!missing(maxProcessors) && maxProcessors==1) {
        out = power_disc(pnull = pnull,
                         rnull = rnull, 
                         vals = vals,
                         ralt = ralt, 
                         param_alt = param_alt, 
                         phat = phat,
                         TS = TS,
                         typeTS = typeTS,
                         TSextra = TSextra,
                         B = B,
                         alpha = alpha)
        rownames(out) = param_alt 
        colnames(out) = methods
        out = cbind(out, chiout)
        if(is.matrix(out) && length(param_alt)==1) out=out[1, ]
        return(round(out,4))
     }  
     if(missing(maxProcessors)) m=parallel::detectCores()-1
     else m=maxProcessors
     cl = parallel::makeCluster(m)
     z=parallel::clusterCall(cl, 
                    power_disc, 
                         pnull = pnull,
                         rnull = rnull, 
                         vals = vals,
                         ralt = ralt, 
                         param_alt = param_alt, 
                         phat = phat,   
                         TS = TS,
                         typeTS = typeTS,
                         TSextra = TSextra,
                         B = c(round(B[1]/m), B[2]),
                         alpha = alpha
                    )
     parallel::stopCluster(cl)
#     Average power of cores    
     out=0*z[[1]]
     for(i in 1:m) out=out+z[[i]]
     out = out/m
     rownames(out) = param_alt 
     colnames(out) = methods
     out = cbind(out, chiout)
     if(is.matrix(out) && length(param_alt)==1) out=out[1, ]
     return(out)
  }
# No parameter estimation    
# critical values of null distributions: 
    if(is.function(phat)) param=phat
    else param=0
    res_rnull=formals(rnull)
    TS_data = matrix(0, B[2], nummethods)    
    for(i in 1:B[2]) {
      x = rnull()
      if(typeTS<2) TS_data[i, ]=TS(x, pnull, phat(x), vals)
      if(typeTS==2) TS_data[i, ]=TS(x, pnull, phat(x), vals, TSextra)
    }  
    crit = apply(TS_data, 2, stats::quantile, probs=1-alpha, na.rm=TRUE)
# power calculations:  
    npar_alt=length(param_alt)
    TS_sim = as.list(1:npar_alt)
    for(i in 1:npar_alt) TS_sim[[i]] = matrix(0, B[1], nummethods)
    for(j in 1:npar_alt) {
      x = ralt(param_alt[j])
      for(i in 1:B[1]) {
        x = ralt(param_alt[j])
        if(typeTS<2) TS_sim[[j]][i, ]=TS(x, pnull, phat(x), vals)
        if(typeTS==2) TS_sim[[j]][i, ]=TS(x, pnull, phat(x), vals, TSextra)
      }
    }
    out = matrix(0, npar_alt, nummethods)
    colnames(out) = methods
    rownames(out) = param_alt
    for(i in seq_along(param_alt)) {
      for(j in 1:nummethods) {
        out[i, j] = sum(TS_sim[[i]][, j]>crit[j])/B[1]
      } 
    }  
    out = cbind(out, chiout)
    if(npar_alt==1) out=out[1, , drop=FALSE]
    out
}

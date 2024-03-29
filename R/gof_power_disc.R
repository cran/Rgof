#' Find the power of various gof tests for discrete data.
#' @param  pnull cumulative distribution function under the null hypothesis
#' @param  rnull  a function to generate data under  null hypothesis
#' @param  vals values of discrete rv.
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat  function to estimate parameters from the data
#' @param  TS user supplied function to find test statistics
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B =c(1000, 1000), number of simulation runs to find power and null distribution
#' @param  nbins =c(100,10), number of bins for chisquare tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  maxProcessors maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  minexpcount =2 minimal expected bin count required
#' @return A numeric matrix of power values.
#' @export 
#' @examples 
#' # Power of tests when null hypothesis specifies a binomial N=10, p=0.5 distribution but 
#' # true data comes from a binomial distribution with success probability 0.55 or 0.6.
#' vals=0:10
#' pnull = function() pbinom(0:10, 10, 0.5)
#' rnull = function() table(c(0:10, rbinom(1000, 10, 0.5)))-1
#' ralt = function(p) table(c(0:10, rbinom(1000, 10, p)))-1
#' gof_power_disc(pnull, rnull, vals, ralt, c(0.515, 0.53), B=c(500, 500))
#' # Power of tests when null hypothesis specifies a binomial N=10 distribution and
#' # p is estimated from the data. 
#' pnull=function(p=0.5)  pbinom(0:10, 10, p)
#' rnull = function(p=0.5) table(c(0:10, rbinom(1000, 10, p)))-1
#' ralt = function(p=0.5) table(c(0:10, rbinom(1000, 10, p)))-1
#' phat = function(x) mean(rep(0:10, x))/10
#' gof_power_disc(pnull, rnull, vals, ralt, phat, param_alt=0.6, B=c(100, 100), maxProcessors = 2)

gof_power_disc=function(pnull, rnull, vals, ralt, param_alt, phat, TS, 
        alpha=0.05, B=c(1000, 1000),  nbins=c(100,10), 
        rate=0, maxProcessors, minexpcount=2.0) {

    if(!missing(phat)) check.functions(pnull, rnull, vals=vals, phat=phat)
    if(missing(phat)) check.functions(pnull, rnull, vals=vals)
    NewTS=FALSE
    if(missing(TS)) TS =  TS_disc
    else {
      if(substr(deparse(TS)[2], 1, 5)==".Call") {
        message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
        maxProcessors=1
      }  
      NewTS=TRUE
    }
    x = ralt(param_alt[1])
    TS_data = TS(x, (1:length(x))/length(x), nm_calc(sum(x)), vals)
    if(NewTS) {
      if(is.null(names(TS_data))) {
        message("result of TS has to be a named vector")
        return(NULL)
      }   
      NamesNewTS=names(TS_data)
    }
    nummethods = length(TS_data)
    methods = c(names(TS_data), "chi large Pearson", "chi small Pearson", 
                paste0(c("chi large LR", "chi small LR"),"-", ifelse(rate==0, "m", "p")))
    if(!missing(phat) && is.function(phat)) {
      if(!missing(maxProcessors) && maxProcessors==1) {
        out = power_disc(pnull = pnull,
                         rnull = rnull, 
                         vals = vals,
                         ralt = ralt, 
                         param_alt = param_alt, 
                         phat = phat,
                         TS = TS,
                         nbins = nbins,
                         rate = rate,
                         B = B,
                         alpha = alpha,
                         minexpcount=minexpcount)
        rownames(out) = param_alt 
        colnames(out) = methods
        if(NewTS) out=out[ , NamesNewTS, drop=FALSE]
        if(is.matrix(out) && length(param_alt)==1) out=out[1, ]
        return(out)
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
                         nbins = nbins,
                         rate = rate,
                         B = c(round(B[1])/m, B[2]),
                         alpha = alpha,
                         minexpcount=minexpcount
                    )
      parallel::stopCluster(cl)
      # Average power of cores    
      out=0*z[[1]]
      for(i in 1:m) out=out+z[[i]]
      rownames(out) = param_alt 
      colnames(out) = methods
      if(NewTS) out=out[ , NamesNewTS, drop=FALSE]
      if(is.matrix(out) && length(param_alt)==1) out=out[1, ]
      return(out/m)
    }
# No parameter estimation    
# critical values of null distributions: 
    if(!missing(phat)) param=phat
    else param=0
    res_pnull=formals(pnull)
    if(length(res_pnull)==0) p = pnull()    
    else p = pnull(param)    
    res_rnull=formals(rnull)
    if(length(res_rnull)==0) nm = nm_calc(sum(rnull()))    
    else nm = nm_calc(sum(rnull(param)))    
    TS_data = matrix(0, B[2], nummethods)    
    for(i in 1:B[2]) {
      x = rnull()
      if(rate>0) nm = nm_calc(x)    
      TS_data[i, ]=TS(x, p, nm, vals)
    }  
    crit = apply(TS_data, 2, stats::quantile, probs=1-alpha)
# power calculations:  
    npar_alt=length(param_alt)
    A = matrix(0, B[1], nummethods)
    TS_sim = list(1:npar_alt)
    for(i in 1:npar_alt) TS_sim[[i]] = A
    A = matrix(0, B[1], 4)
    chi.p.value = list(1:npar_alt)
    for(i in 1:npar_alt) chi.p.value[[i]] = A
    for(i in 1:B[1]) {
      for(j in 1:npar_alt) {
        x = ralt(param_alt[j])
        if(rate>0) nm = nm_calc(sum(x))    
        TS_sim[[j]][i, ] = TS(x, p, nm, vals)      
        chi.p.value[[j]][i, 1:2] = chi_test_disc(x, pnull, param, 
                        nbins, "Pearson", rate, 0, minexpcount)[, 3]
        chi.p.value[[j]][i, 3:4] = chi_test_disc(x, pnull, param, 
                        nbins, "LR", rate, 0, minexpcount)[, 3]
      }  
    }
    out = matrix(0, npar_alt, length(methods))
    colnames(out) = methods
    rownames(out) = param_alt
    for(i in seq_along(param_alt)) {
      for(j in 1:nummethods) {
        out[i, j] = sum(TS_sim[[i]][, j]>crit[j])/B[1]
      } 
      for(j in 1:4) {
        out[i, nummethods+j] = sum(chi.p.value[[i]][, j]<alpha)/B[1]
      }
    }  
    if(NewTS) out=out[ , NamesNewTS, drop=FALSE]    
    if(npar_alt==1) out=out[1, ]
    out
}

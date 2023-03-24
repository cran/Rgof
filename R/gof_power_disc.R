#' Find the power of various gof tests for discrete data.
#' @param  pnull cumulative distribution function under the null hypothesis
#' @param  rnull  a function to generate data under  null hypothesis
#' @param  vals values of discrete rv.
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat  function to estimate parameters from the data
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B =c(1000, 1000), number of simulation runs to find power and null distribution
#' @param  nbins =c(100,10), number of bins for chisquare tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  maxProcessors maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
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

gof_power_disc=function(pnull, rnull, vals, ralt, param_alt, phat, 
        alpha=0.05, B=c(1000, 1000),  nbins=c(100,10), 
        rate=0, maxProcessors) {
    doMethodTS = c("KS", "K","AD", "CvM",  "W", "ZA", "ZK",  "ZC", "Wassp1")    
    doMethodchi = c("chi large Pearson", "chi small Pearson", 
                  paste0(c("chi large LR", "chi small LR"),"-", ifelse(rate==0, "m", "p"))) 
    doMethod = c(doMethodTS, doMethodchi)
    if(!missing(phat)) check.functions(pnull, rnull, vals=vals, phat=phat)
    if(missing(phat)) check.functions(pnull, rnull, vals=vals)
    if(!missing(phat) && is.function(phat)) {
      if(!missing(maxProcessors) && maxProcessors==1) {
        out = power_disc(pnull = pnull,
                         phat = phat,   
                         rnull = rnull, 
                         ralt = ralt, 
                         param_alt = param_alt, 
                         vals = vals,
                         nbins = nbins,
                         rate = rate,
                         B = B,
                         alpha = alpha)
        rownames(out) = param_alt 
        colnames(out) = doMethod
        out = out[, doMethod!="ZK"]  
        return(out)
      }  
      if(missing(maxProcessors)) m=parallel::detectCores()-1
      else m=maxProcessors
      cl = parallel::makeCluster(m)
      z=parallel::clusterCall(cl, 
                    power_disc, 
                         pnull = pnull,
                         phat = phat,   
                         rnull = rnull, 
                         ralt = ralt, 
                         param_alt = param_alt,
                         vals = vals,
                         nbins = nbins,
                         rate = rate,
                         B = c(round(B[1])/m, B[2]),
                         alpha = alpha)
      parallel::stopCluster(cl)
      # Average power of cores    
      out=0*z[[1]]
      for(i in 1:m) out=out+z[[i]]
      rownames(out) = param_alt 
      colnames(out) = doMethod
      out = out[, doMethod!="ZK"]
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
    TS = matrix(0, B[2], 9)    
    for(i in 1:B[2]) {
      x = rnull()
      if(rate>0) nm = nm_calc(x)    
      TS[i, ]=TS_disc(x, p, nm, vals, doMethod)
    }  
    crit = apply(TS, 2, stats::quantile, probs=1-alpha)
# power calculations:  
    npar_alt=length(param_alt)
    A = matrix(0, B[1], 9)
    TS = list(1:npar_alt)
    for(i in 1:npar_alt) TS[[i]] = A
    A = matrix(0, B[1], 4)
    chi.p.value = list(1:npar_alt)
    for(i in 1:npar_alt) chi.p.value[[i]] = A
    for(i in 1:B[1]) {
      for(j in 1:npar_alt) {
        x = ralt(param_alt[j])
        if(rate>0) nm = nm_calc(sum(x))    
        TS[[j]][i, ] = TS_disc(x, p, nm, vals)      
        chi.p.value[[j]][i, 1:2] = chi_test_disc(x, pnull, param, nbins, "Pearson", rate, 0)[, 3]
        chi.p.value[[j]][i, 3:4] = chi_test_disc(x, pnull, param, nbins, "LR", rate, 0)[, 3]
      }  
    }
    pwr = matrix(0, npar_alt, length(doMethod))
    colnames(pwr) = doMethod
    rownames(pwr) = param_alt
    for(i in seq_along(param_alt)) {
      for(j in 1:9) {
        pwr[i, j] = sum(TS[[i]][, j]>crit[j])/B[1]
      } 
      for(j in 1:4) {
        pwr[i, 9+j] = sum(chi.p.value[[i]][, j]<alpha)/B[1]
      }
    }  
    pwr = pwr[, doMethod!="ZK"]
    pwr
}

#' Find the power of various gof tests for continuous data.
#' @param  pnull function to find cdf under  null hypothesis
#' @param  rnull function to generate data under  null hypothesis
#' @param  qnull quantile function (inverse cdf). If missing Wasserstein test can not be done.
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat  function to estimate parameters from the data
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any
#' @param  B =c(1000, 1000), number of simulation runs to find power and null distribution
#' @param  nbins =c(100,10), number of bins for chi square tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  maxProcessors maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @return A numeric matrix of power values.
#' @export 
#' @examples
#' # Power of tests when null hypothesis specifies the standard normal distribution but 
#' # true data comes from a normal distribution with mean different from 0.
#' pnull = function(x) pnorm(x)
#' qnull = function(x) qnorm(x)
#' rnull = function()  rnorm(100)
#' ralt = function(mu)  rnorm(50, mu)
#' gof_power_cont(pnull, rnull, qnull, ralt, c(0.25, 0.5), B=c(500, 500))
#' # Power of tests when null hypothesis specifies normal distribution and 
#' # mean and standard deviation are estimated from the data. 
#' # Example is not run because it takes several minutes.
#' # true data comes from a normal distribution with mean different from 0.
#' pnull = function(x, p=c(0, 1)) pnorm(x, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' qnull = function(x, p=c(0, 1)) qnorm(x, p[1],  ifelse(p[2]>0.001, p[2], 0.001))
#' rnull = function(p=c(0, 1))  rnorm(1000, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' phat = function(x) c(mean(x), sd(x))
#' \donttest{gof_power_cont(pnull, rnull, qnull, ralt, c(0, 1), phat)}

gof_power_cont=function(pnull, rnull, qnull, ralt, param_alt, phat,
        alpha=0.05, Range  =c(-Inf, Inf), B=c(1000, 1000),nbins=c(100,10), 
        rate=0, maxProcessors) {

    doMethod = c("KS", "K","AD", "CvM",  "W", "ZA", "ZK",  "ZC", "Wassp1",    
            "EP large Pearson", "ES large Pearson", "EP small Pearson", "ES small Pearson",
        paste0(c("EP large LR", "ES large LR", "EP small LR", "ES small LR"),"-", ifelse(rate==0, "m", "p")) 
    )
    if(missing(qnull) && missing(phat)) check.functions(pnull, rnull)
    if(!missing(qnull) && missing(phat)) check.functions(pnull, rnull, qnull)
    if(missing(qnull) && !missing(phat)) check.functions(pnull, rnull, phat=phat)
    if(!missing(qnull) && !missing(phat)) check.functions(pnull, rnull, qnull, phat)
    if(missing(qnull)) qnull=function(x, p)  -99
    if(is.infinite(Range[1])) Range[1]=-99999
    if(is.infinite(Range[2])) Range[2]=99999    
    if(!missing(phat) && is.function(phat)) {
# With parameter estimation    
      if(!missing(maxProcessors) && maxProcessors==1) {
         out = power_cont(pnull = pnull,
                          phat = phat,
                          rnull = rnull, 
                          ralt = ralt, 
                          param_alt = param_alt, 
                          qnull = qnull, 
                          nbins = nbins,
                          rate = rate,
                          Range = Range,
                          B = B,
                          alpha = alpha)
         rownames(out) = param_alt 
         colnames(out) = doMethod  
         return(out)
    }  
    if(missing(maxProcessors)) m=parallel::detectCores()-1
    else m=maxProcessors
    cl = parallel::makeCluster(m)
    z=parallel::clusterCall(cl, 
                          power_cont, 
                              pnull = pnull,
                              phat = phat, 
                              rnull = rnull, 
                              ralt = ralt, 
                              param_alt = param_alt, 
                              qnull = qnull,
                              nbins = nbins,
                              rate = rate,
                              Range = Range,
                              B = c(round(B[1])/m, B[2]),
                              alpha = alpha)
      parallel::stopCluster(cl)
      # Average power of cores
      out=0*z[[1]]
      for(i in 1:m) out=out+z[[i]]
      colnames(out)=doMethod
      rownames(out)=param_alt 
      return(out/m)
    }
# No parameter estimation   

# critical values of null distributions: 
    if(!missing(phat)) param=phat
    else param=0
    res_pnull=formals(pnull)
    res_rnull=formals(rnull)
    TS = matrix(0, B[2], 9)    
    for(i in 1:B[2]) {
      if(length(res_rnull)==0) x=rnull()
      else x=rnull(x, param)
      if(length(res_pnull)==1) Fx=pnull(x)
      else Fx=pnull(x, param)
      TS[i, ]=TS_cont(x, Fx, 0, qnull, doMethod)
    }  
    crit = apply(TS, 2, stats::quantile, probs=1-alpha)
# power calculations:
    npar_alt = length(param_alt)
    A = matrix(0, B[1], 9)
    TS = list(1:npar_alt)
    for(i in 1:npar_alt) TS[[i]] = A
    A = matrix(0, B[1], 8)
    chi.p.value = list(1:npar_alt)
    Range=ifelse(is.infinite(Range),-99, Range)
    for(i in 1:npar_alt) chi.p.value[[i]] = A
    for(i in 1:B[1]) {
        for(j in 1:npar_alt) {
           x = ralt(param_alt[j])
           if(length(res_pnull)==1) Fx=pnull(x)
           else Fx=pnull(x, param)
           TS[[j]][i, ] = TS_cont(x, Fx, param, qnull, doMethod)
           chi.p.value[[j]][i, 1:4] = 
              chi_test_cont(x, pnull, param, "Pearson", rate, nbins, Range, 0)[, 3]
           chi.p.value[[j]][i, 5:8] = 
             chi_test_cont(x, pnull, param, "LR", rate, nbins, Range, 0)[, 3]           
        }  
    }

    pwr = matrix(0, npar_alt, 17)
    colnames(pwr) = doMethod
    rownames(pwr) = param_alt
    for(i in 1:npar_alt) {
       for(j in 1:9) {
          pwr[i, j] = sum(TS[[i]][, j]>crit[j])/B[1]
       } 
       for(j in 1:8) {
          pwr[i, 9+j] = sum(chi.p.value[[i]][,j]<alpha)/B[1]
      }         
    }

    pwr
}

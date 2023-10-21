#' Find the power of various gof tests for continuous data.
#' @param  pnull function to find cdf under  null hypothesis
#' @param  rnull function to generate data under  null hypothesis
#' @param  qnull quantile function (inverse cdf). If missing Wasserstein test can not be done.
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat  function to estimate parameters from the data
#' @param  TS user supplied function to find test statistics
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any
#' @param  B =c(1000, 1000), number of simulation runs to find power and null distribution
#' @param  nbins =c(100,10), number of bins for chi square tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  maxProcessors maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  minexpcount =2 minimal expected bin count required
#' @return A numeric matrix of power values.
#' @export 
#' @examples
#' # Power of tests when null hypothesis specifies the standard normal distribution but 
#' # true data comes from a normal distribution with mean different from 0.
#' pnull = function(x) pnorm(x)
#' qnull = function(x) qnorm(x)
#' rnull = function()  rnorm(50)
#' ralt = function(mu)  rnorm(50, mu)
#' gof_power_cont(pnull, rnull, qnull, ralt, c(0.25, 0.5), B=c(500, 500))
#' # Power of tests when null hypothesis specifies normal distribution and 
#' # mean and standard deviation are estimated from the data. 
#' # Example is not run because it takes several minutes.
#' # true data comes from a normal distribution with mean different from 0.
#' pnull = function(x, p=c(0, 1)) pnorm(x, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' qnull = function(x, p=c(0, 1)) qnorm(x, p[1],  ifelse(p[2]>0.001, p[2], 0.001))
#' rnull = function(p=c(0, 1))  rnorm(50, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' phat = function(x) c(mean(x), sd(x))
#' \donttest{gof_power_cont(pnull, rnull, qnull, ralt, c(0, 1), phat, B=c(200, 200), maxProcessor=2)}

gof_power_cont=function(pnull, rnull, qnull, ralt, param_alt, phat, TS, 
        alpha=0.05, Range  =c(-Inf, Inf), B=c(1000, 1000),nbins=c(100,10), 
        rate=0, maxProcessors, minexpcount=2.0) {
  # Do some sanity checks and setups 
  if(missing(qnull) && missing(phat)) check.functions(pnull, rnull)
  if(!missing(qnull) && missing(phat)) check.functions(pnull, rnull, qnull)
  if(missing(qnull) && !missing(phat)) check.functions(pnull, rnull, phat=phat)
  if(!missing(qnull) && !missing(phat)) check.functions(pnull, rnull, qnull, phat)
  if(missing(qnull)) qnull=function(x, p)  -99
  if(is.infinite(Range[1])) Range[1]=-99999
  if(is.infinite(Range[2])) Range[2]=99999  
  NewTS=FALSE
  if(missing(TS)) TS =  TS_cont
  else { # can't do parallel processing
    if(substr(deparse(TS)[2], 1, 5)==".Call") {
      message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
      maxProcessors=1
    }  
    NewTS=TRUE
  }
  # Find out how many (non-chi-square) methods are implemented and what they are called
  x=ralt(param_alt[1])
  TS_data=TS(x, (1:length(x))/(length(x)+1), 1, qnull)
  if(NewTS) {
     if(is.null(names(TS_data))) {
        message("result of TS has to be a named vector")
        return(NULL)
     }   
     NamesNewTS=names(TS_data)
  }  
  nummethods = length(TS_data)
  methods = c(names(TS_data),    
            "EP large Pearson", "ES large Pearson", "EP small Pearson", "ES small Pearson",
        paste0(c("EP large LR", "ES large LR", "EP small LR", "ES small LR"),"-", ifelse(rate==0, "m", "p")) 
    )

    if(!missing(phat) && is.function(phat)) {
      
# With parameter estimation    
      if(!missing(maxProcessors) && maxProcessors==1) {
         out = power_cont(pnull = pnull,
                          rnull = rnull,
                          qnull = qnull,  
                          ralt = ralt, 
                          param_alt = param_alt, 
                          phat = phat,
                          TS = TS,
                          nbins = nbins,
                          rate = rate,
                          Range = Range,
                          B = B,
                          alpha = alpha,
                          minexpcount=minexpcount
                          )
         rownames(out) = param_alt 
         colnames(out) = methods
         if(NewTS) out=out[ , NamesNewTS, drop=FALSE]
         if(is.matrix(out) & nrow(out)==1) out=out[1, ]
         return(out)
    }  
    if(missing(maxProcessors)) m=parallel::detectCores()-1
    else m=maxProcessors
    cl = parallel::makeCluster(m)
    z=parallel::clusterCall(cl, 
                          power_cont, 
                              pnull = pnull,
                              rnull = rnull,
                              qnull = qnull,  
                              ralt = ralt, 
                              param_alt = param_alt, 
                              phat = phat,
                              TS = TS,
                              nbins = nbins,
                              rate = rate,
                              Range = Range,
                              B = c(round(B[1])/m, B[2]),
                              alpha = alpha,
                              minexpcount=minexpcount
                          )
      parallel::stopCluster(cl)
      # Average power of cores
      out=0*z[[1]]
      for(i in 1:m) out=out+z[[i]]
      colnames(out)=methods
      rownames(out)=param_alt 
      if(NewTS) out=out[ , NamesNewTS, drop=FALSE]
      if(is.matrix(out) & nrow(out)==1) out=out[1, ]
      return(out/m)
    }
# No parameter estimation   

# critical values of null distributions: 
    if(!missing(phat)) param=phat
    else param=0
    res_pnull=formals(pnull)
    res_rnull=formals(rnull)
    TS_data = matrix(0, B[2], nummethods)    
    for(i in 1:B[2]) {
      if(length(res_rnull)==0) x=rnull()
      else x=rnull(x, param)
      if(length(res_pnull)==1) Fx=pnull(x)
      else Fx=pnull(x, param)
      TS_data[i, ]=TS(x, Fx, 0, qnull)
    }  
    crit = apply(TS_data, 2, stats::quantile, probs=1-alpha)
# power calculations:
    npar_alt = length(param_alt)
    A = matrix(0, B[1], nummethods)
    TS_alt = list(1:npar_alt)
    for(i in 1:npar_alt) TS_alt[[i]] = A
    A = matrix(0, B[1], 8)
    chi.p.value = list(1:npar_alt)
    Range=ifelse(is.infinite(Range),-99, Range)
    for(i in 1:npar_alt) chi.p.value[[i]] = A
    for(i in 1:B[1]) {
        for(j in 1:npar_alt) {
           x = ralt(param_alt[j])
           if(length(res_pnull)==1) Fx=pnull(x)
           else Fx=pnull(x, param)
           TS_alt[[j]][i, ] = TS(x, Fx, param, qnull)
           chi.p.value[[j]][i, 1:4] = 
              chi_test_cont(x, pnull, param, "Pearson", rate, nbins, Range, 0, minexpcount)[, 3]
           chi.p.value[[j]][i, 5:8] = 
             chi_test_cont(x, pnull, param, "LR", rate, nbins, Range, 0, minexpcount)[, 3]           
        }  
    }

    out = matrix(0, npar_alt, nummethods+8)
    colnames(out) = methods
    rownames(out) = param_alt
    for(i in 1:npar_alt) {
       for(j in 1:nummethods) {
          out[i, j] = sum(TS_alt[[i]][, j]>crit[j])/B[1]
       } 
       for(j in 1:8) {
          out[i, nummethods+j] = sum(chi.p.value[[i]][,j]<alpha)/B[1]
      }         
    }
    if(NewTS) out=out[ , NamesNewTS, drop=FALSE]
    if(npar_alt==1) out=out[1, ]
    out
}

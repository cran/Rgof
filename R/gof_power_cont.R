#' Find the power of various gof tests for continuous data.
#' @param  pnull function to find cdf under  null hypothesis
#' @param  rnull function to generate data under  null hypothesis
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  w (Optional) function to calculate weights, returns -99 if no weights
#' @param  phat =function(x) -99, function to estimate parameters from the data, or -99 if no parameters aare estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  TSextra =NA, list provided to TS
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any
#' @param  B =c(1000, 1000), number of simulation runs to find power and null distribution
#' @param  nbins =c(100,10), number of bins for chi square tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  maxProcessors maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat =TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @return A numeric matrix of power values.

gof_power_cont=function(pnull, rnull, ralt, param_alt, 
        w=function(x) -99, phat=function(x) -99, TS, TSextra=NA, 
        alpha=0.05, Range  =c(-Inf, Inf), B=c(1000, 1000),nbins=c(100,10), 
        rate=0, maxProcessors, minexpcount=5.0, ChiUsePhat=TRUE) {
  
  x = ralt(param_alt[1])  
  WithWeights = TRUE
  if(length(formals(w))==1 && w(x[1])==-99) WithWeights = FALSE
  WithEstimation = TRUE
  if(phat(x)[1]==-99) WithEstimation = FALSE

  if(any(is.na(TSextra))) TSextra = list(p=phat(x))
  else TSextra = c(TSextra, p=phat)
  Noqnull = FALSE
  if( !("qnull" %in% names(TSextra)) ) {
    Noqnull = TRUE
    qnull=function(x) abs(x)/sum(abs(x))
    TSextra = c(TSextra, qnull=qnull)
  }  
  else qnull = TSextra$qnull
  if(missing(TS)) {
    if(!WithWeights) { #data is not weighted
      typeTS=1
      TS = TS_cont
      TS_data = TS(x, pnull, phat(x), function(x) abs(x)/max(x))
    }
    else {
      typeTS=2
      TS = TSw_cont
      if(length(formals(w))==1) wx=w(x)
      else wx=w(x, phat(x))
      TS_data = TS(x, pnull, phat(x), wx)
      doMethods = names(TS_data)
    }
  }   
  else {
    # can't do parallel processing if TS written in C/C++
    if(substr(deparse(TS)[2], 1, 5)==".Call") {
      message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
      maxProcessors=1
    }
    if(length(formals(TS))==3) {
      typeTS=3
      TS_data = TS(x, pnull, phat(x))
    }
    if(length(formals(TS))==4) {
      typeTS=4
      TS_data = TS(x, pnull, phat(x), TSextra)
    }
    if(length(formals(TS))>4) {
      message("TS should have either 3 or 4 arguments")
      return(NULL)
    }
    if(is.null(names(TS_data))) {
      message("result of TS has to be a named vector")
      return(NULL)
    }
  }
  nummethods = length(TS_data)
  methods = names(TS_data)
# Do chi square tests if built-in TS is used. Don't run chi square if weights are present.  
  chiout=NULL
  if(typeTS<=2 && !WithWeights) { #Run chi square tests
    if(is.infinite(Range[1])) Range[1]=-99999
    if(is.infinite(Range[2])) Range[2]=99999  
    chiout = chi_power_cont(pnull=pnull, 
                            ralt = ralt, 
                            param_alt = param_alt,                             
                            qnull = ifelse(Noqnull, NA, qnull), 
                            phat = phat, 
                            w = w,
                            alpha = alpha, 
                            Range = Range, 
                            B= B[1], 
                            nbins = nbins, 
                            rate = rate, 
                            minexpcount = minexpcount,
                            ChiUsePhat=ChiUsePhat)  
    if(WithWeights) chiout = chiout[, 1:4]
  }
# Now the other tests
  if(WithEstimation) {
# With parameter estimation 
    if(missing(maxProcessors)) m=parallel::detectCores()-1
    else m=maxProcessors
    if(m==1) {
         out = power_cont(pnull = pnull,
                          rnull = rnull,
                          qnull = qnull,  
                          ralt = ralt, 
                          param_alt = param_alt, 
                          phat = phat,
                          w = w,
                          TS = TS,
                          typeTS = typeTS,
                          TSextra = TSextra,
                          B = B,
                          alpha = alpha
                          )
         rownames(out) = param_alt 
         colnames(out) = methods
         if(Noqnull) out = out[ ,colnames(out)!="Wassp1", drop=FALSE] 
         out = cbind(out, chiout)
         if(is.matrix(out) & nrow(out)==1) out=out[1, ]
         return(out)
    }  
    cl = parallel::makeCluster(m)
    z=parallel::clusterCall(cl, 
                    power_cont, 
                          pnull = pnull,
                          rnull = rnull,
                          qnull = qnull,  
                          ralt = ralt, 
                          param_alt = param_alt, 
                          w = w,
                          phat = phat,
                          TS = TS,
                          typeTS = typeTS,
                          TSextra = TSextra,
                          B = c(round(B[1]/m), B[2]),
                          alpha = alpha
                    )
      parallel::stopCluster(cl)
      # Average power of cores
      out=0*z[[1]]
      for(i in 1:m) out=out+z[[i]]
      colnames(out)=methods
      rownames(out)=param_alt 
      if(Noqnull) out = out[ ,colnames(out)!="Wassp1", drop=FALSE]
      out = cbind(out/m, chiout)
      if(is.matrix(out) & nrow(out)==1) out=out[1, ]
      return(out)
  }
  
# No parameter estimation   

# critical values of null distributions: 
    if(!is.function(phat)) param=phat
    else param=0
    res_pnull=formals(pnull)
    res_rnull=formals(rnull)
    TS_data = matrix(0, B[2], nummethods)    
    for(i in 1:B[2]) {
      if(length(res_rnull)==0) x=rnull()
      else x=rnull(x, param)
      if(typeTS==1) TS_data[i, ]=TS(x, pnull, phat(x), qnull)
      if(typeTS==2) TS_data[i, ]=TS(x, pnull, phat(x), w(x))
      if(typeTS==3) TS_data[i, ]=TS(x, pnull, phat(x))
      if(typeTS==4) TS_data[i, ]=TS(x, pnull, phat(x), TSextra)
    }  
    crit = apply(TS_data, 2, stats::quantile, probs=1-alpha)
# power calculations:
    npar_alt = length(param_alt)
    A = matrix(0, B[1], nummethods)
    TS_alt = list(1:npar_alt)
    for(i in 1:npar_alt) TS_alt[[i]] = A
    for(i in 1:B[1]) {
        for(j in 1:npar_alt) {
           x = ralt(param_alt[j])
           if(typeTS==1) TS_alt[[j]][i, ]=TS(x, pnull, phat(x), qnull)
           if(typeTS==2) TS_alt[[j]][i, ]=TS(x, pnull, phat(x), w(x))
           if(typeTS==3) TS_alt[[j]][i, ]=TS(x, pnull, phat(x))
           if(typeTS==4) TS_alt[[j]][i, ]=TS(x, pnull, phat(x), TSextra)
        }  
    }

    out = matrix(0, npar_alt, nummethods)
    colnames(out) = methods
    rownames(out) = param_alt
    for(i in 1:npar_alt) {
       for(j in 1:nummethods) {
          out[i, j] = sum(TS_alt[[i]][, j]>crit[j])/B[1]
       } 
    }
    out = cbind(out, chiout)  
    if(Noqnull) out = out[ ,colnames(out)!="Wassp1", drop=FALSE]
    if(npar_alt==1) out=out[1, , drop=FALSE]
    out
}

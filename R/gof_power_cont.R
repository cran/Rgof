#' Find the power of various gof tests for continuous data.
#' @param  pnull function to find cdf under  null hypothesis
#' @param  rnull function to generate data under  null hypothesis
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  w =function(x) -99, function to calculate weights, returns -99 if no weights
#' @param  phat =function(x) -99, function to estimate parameters from the data, or -99 if no parameters aare estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  TSextra =NA, list provided to TS
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any
#' @param  B =1000 number of simulation runs
#' @param  nbins =c(100,10), number of bins for chi square tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  maxProcessor maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat =TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @return A numeric matrix of power values.

gof_power_cont=function(pnull, rnull, ralt, param_alt, 
        w=function(x) -99, phat=function(x) -99, TS, TSextra=NA, 
        alpha=0.05, Range  =c(-Inf, Inf), B=1000,nbins=c(100,10), 
        rate=0, maxProcessor, minexpcount=5.0, ChiUsePhat=TRUE) {
  
  x = ralt(param_alt[1])  
  WithWeights = TRUE
  if(length(formals(w))==1 && w(x[1])==-99) WithWeights = FALSE
  WithEstimation = TRUE
  if(phat(x)[1]==-99) WithEstimation = FALSE
  if(any(is.na(TSextra))) TSextra = list(p=phat(x))
  else TSextra = c(TSextra, p=phat)
  TSextra = c(TSextra, w=w)
  Noqnull = FALSE
  if( !("qnull" %in% names(TSextra)) ) {
    Noqnull = TRUE
    qnull=function(x) -99
    TSextra = c(TSextra, qnull=qnull)
  }  
  else qnull = TSextra$qnull
  if(missing(TS)) {
    if(!WithWeights) { #data is not weighted
      typeTS=1
      TS = TS_cont
      TS_data = TS(x, pnull, phat(x), qnull)
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
      maxProcessor=1
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
                            B= B, 
                            nbins = nbins, 
                            rate = rate, 
                            minexpcount = minexpcount,
                            ChiUsePhat=ChiUsePhat)  
    if(WithWeights) chiout = chiout[, 1:4]
  }
# Now the other tests
  out = power_cont_R(pnull = pnull,
                    rnull = rnull,  
                    ralt = ralt, 
                    param_alt = param_alt,
                    phat = phat,
                    TS = TS,
                    typeTS = typeTS,
                    TSextra = TSextra,
                    alpha = alpha,
                    B = B,
                    maxProcessor =  maxProcessor
        )
   rownames(out) = param_alt 
   colnames(out) = methods
   if(Noqnull) out = out[ ,colnames(out)!="Wassp1", drop=FALSE] 
   out = cbind(out, chiout)
   if(is.matrix(out) & nrow(out)==1) out=out[1, ]
   out
}

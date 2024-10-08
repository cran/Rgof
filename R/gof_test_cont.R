#' This function performs a number of gof tests for continuous data
#' @param  x data set
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  w (Optional) function to calculate weights, returns -99 if no weights
#' @param  phat =function(x) -99, function to estimate parameters from the data, or -99 if no parameters aare estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  TSextra =NA, list passed to TS, if desired
#' @param  nbins =c(50, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any, for chi-square tests
#' @param  B   =5000  number of simulation runs
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat =TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  maxProcessors =1, number of processors to use in parallel processing. If missing single processor is used.
#' @param  doMethods Methods to include in tests
#' @return A list with vectors of test statistics and p.values

gof_test_cont <- function(x, pnull,  rnull, w=function(x) -99, phat=function(x) -99, 
                          TS, TSextra=NA, nbins=c(50, 10), rate=0, 
                          Range=c(-Inf, Inf), B=5000,  minexpcount=5.0, 
                          ChiUsePhat=TRUE, maxProcessors=1, doMethods="all") {
  # Are weights present?
  WithWeights = TRUE
  if(length(formals(w))==1 & w(x[1])==-99) WithWeights = FALSE
  if(any(is.na(TSextra))) TSextra = list(p=phat(x))
  else TSextra = c(TSextra, p=phat)
  Noqnull = FALSE
  if( !("qnull" %in% names(TSextra)) ) {
    Noqnull = TRUE
    qnull=function(x, p=0) rep(-99,length(x))
    TSextra = c(TSextra, qnull=qnull)
  }  
  else qnull = TSextra$qnull
  if(missing(TS)) {
     nn = 1:length(x)/length(x)
     if(!WithWeights) { #data is not weighted
       typeTS=1
       TS = TS_cont
       TS_data = TS(x, nn, 0, function(x) abs(x)/max(x))
     }
     else {
       typeTS=2
       TS = TSw_cont
       TS_data = TS(x, nn, w(x))
       doMethods = names(TS_data)
     }
  }   
  else {
    # can't do parallel processing if TS written in C/C++
    if(substr(deparse(TS)[2], 1, 5)==".Call") {
      message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
      maxProcessors=1
    }
    if(length(formals(TS))==2) {
       typeTS=3
       TS_data = TS(x, (1:length(x))/(length(x)+1))
    }
    if(length(formals(TS))==3) {
      typeTS=4
      TS_data = TS(x, (1:length(x))/(length(x)+1), TSextra)
    }
    if(length(formals(TS))>3) {
      message("TS should have either 2 or 3 arguments")
      return(NULL)
    }
    if(is.null(names(TS_data))) {
      message("result of TS has to be a named vector")
      return(NULL)
    }
  }
  
  if(maxProcessors==1)
      out = gof_cont(x, pnull, rnull, qnull, w, phat, TS, typeTS, TSextra, B)
  else {
      m=maxProcessors
      cl = parallel::makeCluster(m)
      z=parallel::clusterCall(cl, 
                                gof_cont, 
                                x = x,
                                pnull = pnull,
                                rnull = rnull,
                                qnull = qnull,
                                w = w,
                                phat = phat,
                                TS = TS,
                                typeTS = typeTS,
                                TSextra = TSextra,
                                B = B/m
      )
      parallel::stopCluster(cl)
      # Average power of cores
      tmp=0*z[[1]]
        for(i in 1:m) tmp=tmp+z[[i]]
        out = tmp/m  
  }
  if(typeTS>2) return(list(statistics=out[1, ], p.values=out[2, ]) )
  # do chi square tests
  if(is.infinite(Range[1])) Range[1]=-99999
  if(is.infinite(Range[2])) Range[2]=99999
  outchi = t(chi_test_cont(x, pnull, w, phat, 
          ifelse(Noqnull, NA, qnull),
          nbins, rate, Range, minexpcount, ChiUsePhat)[,c(1, 2)])                        
  if(WithWeights) outchi=outchi[, 1:4]
  out = cbind(out, outchi) 
  if(doMethods[1]=="Default")        
     out = out[ ,c("W", "ZK", "ZC", "Wassp1", "EP-s-P","ES-s-P")]
  if(doMethods[1]!="Default" & doMethods[1]!="all" & typeTS==0) 
     out = out[ ,doMethods, drop=FALSE]
  if(Noqnull) out = out[ ,colnames(out)!="Wassp1"]
  out = round(out, 4)
  list(statistics=out[1, ], p.values=out[2, ])
}

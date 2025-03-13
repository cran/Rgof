#' This function performs a number of gof tests for discrete data.
#' @param  x data set (the counts)
#' @param  pnull  cumulative distribution function under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  vals a vector of values of discrete random variables 
#' @param  phat =function(x) -99, function to estimate parameters from the data, or -99 if no parameters aare estimated
#' @param  TS =NA, user supplied function to find test statistics
#' @param  TSextra =NA, list passed to TS, if desired
#' @param  nbins =c(50, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  B   =5000  number of simulation runs
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  maxProcessor =1, number of processors to use in parallel processing. If missing single processor is used.
#' @param  doMethods Methods to include in tests
#' @return A numeric matrix of test statistics and p.values

gof_test_disc <- function(x, pnull, rnull, vals, phat=function(x) -99, 
                          TS, TSextra=NA,  nbins=c(50, 10), rate=0, 
                          B=5000, minexpcount=5.0, ChiUsePhat=TRUE,
                          maxProcessor=1, doMethods="Default") {

  if(any(is.na(TSextra))) TSextra = list(p=phat(x))
  else TSextra = c(TSextra, p=phat)
  if(missing(TS)) { # use built-in tests
    typeTS = 5
    TS = TS_disc
    TS_data = TS(x, pnull, phat(x), vals)
  }  
  else {
    # can't do parallel processing if TS written in C/C++
    if(substr(deparse(TS)[2], 1, 5)==".Call") {
      message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
      maxProcessor=1
    }
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
  }
  nummethods = length(TS_data)
  out = matrix(0, 2, nummethods)    
  colnames(out) = names(TS_data)

  if(maxProcessor==1)
    out = gof_disc(x, pnull, rnull, vals, phat, TS, typeTS, TSextra, rate, B)
  else {
    m = maxProcessor
    cl = parallel::makeCluster(m)
    z=parallel::clusterCall(cl, 
                            gof_disc, 
                            x = x,
                            pnull = pnull,
                            rnull = rnull,
                            vals = vals,  
                            phat = phat,
                            TS = TS,
                            typeTS = typeTS,
                            TSextra=TSextra,
                            rate=rate,
                            B = B/m
    )
    parallel::stopCluster(cl)
    #  Average power of cores
    tmp=0*z[[1]]
    for(i in 1:m) tmp=tmp+z[[i]]
    out[ ,1:nummethods] = tmp/m  
  }
  if(typeTS>0) return(list(statistics=out[1, ], p.values=out[2, ]))
# do chi square tests
  chiout = t(chi_test_disc(x, pnull, phat, 
              nbins, rate, minexpcount, ChiUsePhat)[,1:2])
  out = cbind(out, chiout)
  if(doMethods[1]=="Default")        
     out = out[ ,c("K", "AD", "ZA", "ZC")]
  if(doMethods[1]!="Default" & doMethods[1]!="all") 
     out = out[ ,doMethods, drop=FALSE]
  out = round(out, 4)
  list(statistics=out[1, ], p.values=out[2, ])
}

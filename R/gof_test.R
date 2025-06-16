#' Tests for the univariate goodness-of-fit problem
#' 
#' This function runs a number of goodness-of-fit tests using Rcpp and parallel computing.
#' 
#' For details on the usage of this routine consult the vignette with vignette("Rgof","Rgof")
#'
#' @param  x data set
#' @param  vals =NA, values of discrete RV, or NA if data is continuous
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  w (Optional) function to calculate weights, returns -99 if no weights
#' @param  phat  =function(x) -99, function to estimate parameters from the data, or -99 if no parameters are estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  TSextra =NA, list passed to TS, if desired, or NA
#' @param  nbins =c(100, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any, for chi-square tests
#' @param  B   =5000  number of simulation runs
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  maxProcessor =1, number of processors to use in parallel processing. 
#' @param  doMethods ="all", a vector of codes for the methods to include or all of them. 
#' @return A list with vectors of test statistics and p.values
#' @export
#' @examples
#' # Tests to see whether data comes from a standard normal distribution.
#' pnull = function(x) pnorm(x)
#' rnull = function()  rnorm(100)
#' x = rnorm(100)
#' gof_test(x, NA, pnull, rnull, B=500)
#' # Tests to see whether data comes from a normal distribution with standard deviation 1 
#' # and the mean estimated.
#' pnull=function(x, m) pnorm(x, m)
#' rnull=function(m) rnorm(100, m)
#' TSextra = list(qnull=function(x, m=0) qnorm(x, m), 
#'           pnull=function(x, m=0) pnorm(x, m), phat=function(x) mean(x))
#' phat=function(x) mean(x)
#' x = rnorm(100, 1, 2)
#' gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra, B=500)
#' # Tests to see whether data comes from a binomial (10, 0.5) distribution.
#' vals=0:10
#' pnull = function() pbinom(0:10, 10, 0.5)
#' rnull = function() table(c(0:10, rbinom(1000, 10, 0.5)))-1
#' x = rnull() 
#' gof_test(x, vals, pnull, rnull, doMethods="all", B=500)
#' # Tests to see whether data comes from a binomial distribution with 
#' # the success probability estimated from the data.
#' pnull = function(p=0.5) pbinom(0:10, 10, ifelse(p>0&&p<1, p, 0.001))
#' rnull = function(p=0.5) table(c(0:10, rbinom(1000, 10, 
#'                   ifelse(p>0&&p<1, p, 0.001))))-1
#' phat=function(x) mean(rep(0:10,x))/10 
#' gof_test(x, vals, pnull, rnull, phat=phat, B=500) 
#'
gof_test <- function(x, vals= NA, pnull, rnull, 
                    w=function(x) -99, phat=function(x) -99, 
                    TS, TSextra=NA, nbins=c(50, 10), rate=0, 
                    Range=c(-Inf, Inf), B=5000,  minexpcount=5.0,  
                    ChiUsePhat=TRUE, maxProcessor, doMethods="all") {

  if(any(is.na(vals))) {
    Continuous=TRUE
    dta=list(x=x)
    check.functions(pnull, rnull, phat, x=x)
  }  
  else {
    Continuous=FALSE
    dta=list(x=x, vals=vals)
    check.functions(pnull, rnull, phat, vals, x)
  }  
  if(any(is.na(TSextra))) TSextra=list(pnull=pnull, phat=phat, 
                      w=w, Continuous=Continuous)
  else TSextra = c(TSextra, pnull=pnull, phat=phat, 
                   w=w, Continuous=Continuous)
  Noqnull=FALSE
  if(!("qnull" %in% names(TSextra))) {
    Noqnull=TRUE
    TSextra=c(TSextra, qnull=function(x) -99)
  }   
  WithWeights = TRUE
  if(length(formals(w))==1) {
    if(w(x[1])==-99) WithWeights = FALSE
  }  
  if(length(x)>10000 && maxProcessor==1)
      message("Consider using parallel processing with maxProcessor= (your number of cores)")
  if(test_methods(doMethods, Continuous, WithWeights)) return(NULL)
# adjust number of bins to account for parameter estimation
   if(abs(phat(x)[1]+99)<0.001) nbins=nbins+length(phat(x)) 
   if(missing(TS)) {
    if(Continuous) {
      if(!WithWeights) { #data is not weighted
        typeTS=1
        TS = TS_cont
      }
      else {
        typeTS=2
        TS = TSw_cont
      }
    }
    else {
      typeTS = 5
      TS = TS_disc
    }
  }   
  else {
    # can't do parallel processing if TS written in C/C++
    if(substr(deparse(TS)[2], 1, 5)==".Call") {
      message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
      maxProcessor=1
    }
    if(Continuous) {
      if(length(formals(TS))>4) {
        message("TS for continuous data should have either 3 or 4 arguments")
        return(NULL)
      }  
      typeTS=length(formals(TS))
    }  
    else {
      if(length(formals(TS))>6) {
        message("TS for discrete datashould have either 4 or 5 arguments")
        return(NULL)
      }
      typeTS=length(formals(TS))+1
    }
  }
  TS_data=calcTS(dta, TS, typeTS, TSextra)  
  if(is.null(names(TS_data))) {
    message("result of TS has to be a named vector")
    return(NULL)
  }
  if(missing(maxProcessor)) 
    maxProcessor=parallel::detectCores(logical = FALSE)-1
  if(maxProcessor>1) {
    tm=timecheck(dta, TS, typeTS, TSextra)
    if(tm*B<20) {
      maxProcessor=1
      message("maxProcessor set to 1 for faster computation")
    }
    else message(paste("Using ",maxProcessor," cores.."))
  }
  if(maxProcessor==1) 
     outTS=gof_test_C(dta, rnull, TS, typeTS, TSextra, B) 
     
  else {
     cl = parallel::makeCluster(maxProcessor)
     z=parallel::clusterCall(cl, gof_test_C, 
             dta, rnull, TS, typeTS, TSextra, 
             B=round(B/maxProcessor)
     )
     parallel::stopCluster(cl)
     #  Average power of cores
     tmp=0*z[[1]]
     for(i in 1:maxProcessor) tmp=tmp+z[[i]]
     outTS = tmp/maxProcessor  
  }
  # do chi square tests
  if(typeTS==1) {
    if(is.infinite(Range[1])) Range[1]=-99999
    if(is.infinite(Range[2])) Range[2]=99999
    outchi = t(chi_test_cont(x, pnull, w, phat, 
                    ifelse(Noqnull, NA, TSextra$qnull),
                    nbins, rate, Range, minexpcount, ChiUsePhat)[,c(1, 2)])
  }
  if(typeTS==5) {
    outchi = t(chi_test_disc(x, pnull, phat, 
                    nbins, rate, minexpcount, ChiUsePhat)[,1:2])
  }
  if(typeTS==1 | typeTS==5) {
    out=list(statistics=c(outTS[1, ], outchi[1, ]), 
             p.values=c(outTS[2, ], outchi[2, ]))
  }           
  else out=list(statistics=outTS[1, ], p.values=outTS[2, ])
  if(doMethods[1]!="all") {
    out[[1]]=out[[1]][doMethods]
    out[[2]]=out[[2]][doMethods]
  }  
  if(Noqnull) {
    out[[1]]=out[[1]][names(out[[1]])!="Wassp1"]
    out[[2]]=out[[2]][names(out[[2]])!="Wassp1"]
  }
  # make output look nice
  signif.digits(out)
}

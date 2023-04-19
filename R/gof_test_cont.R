#' This function performs a number of gof tests for continuous data
#' @param  x data set
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  qnull  routine to calculate quantiles under null hypothesis
#' @param  phat   function to estimate parameters from the data
#' @param  TS user supplied function to find test statistics
#' @param  nbins =c(100, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any, for chi-square tests
#' @param  B   =5000  number of simulation runs
#' @param  doMethod Methods to include in tests
#' @return A list with vectors of test statistics and p values
#' @export
#' @examples
#' # Tests to see whether data comes from a standard normal distribution.
#' pnull = function(x) pnorm(x)
#' qnull = function(x) qnorm(x)
#' rnull = function()  rnorm(100)
#' x = rnorm(100) 
#' gof_test_cont(x, pnull, rnull, qnull, doMethod="all")
#' # Tests to see whether data comes from a normal distribution with 
#' # mean and standard deviation estimated from the data.
#' pnull = function(x, p=c(0, 1)) pnorm(x, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' qnull = function(x, p=c(0, 1)) qnorm(x, p[1],  ifelse(p[2]>0.001, p[2], 0.001))
#' rnull = function(p=c(0, 1))  rnorm(100, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' phat = function(x) c(mean(x), sd(x))
#' gof_test_cont(x, pnull, rnull, qnull, phat)
#' 
gof_test_cont <- function(x, pnull,  rnull, qnull, phat, TS, 
       nbins=c(100, 10), rate=0, Range=c(-Inf, Inf), B=5000,  doMethod="Default") {

  if(missing(qnull) && missing(phat)) check.functions(pnull, rnull)
  if(!missing(qnull) && missing(phat)) check.functions(pnull, rnull, qnull)
  if(missing(qnull) && !missing(phat)) check.functions(pnull, rnull, phat=phat)
  if(!missing(qnull) && !missing(phat)) check.functions(pnull, rnull, qnull, phat)  

  if(missing(TS)) TS = TS_cont
# if phat or qnull is missing, create it
  Minimize = 1
  if(missing(phat)) {
    Minimize = 0
    phat = function(x) 0
  }
  if(missing(qnull)) {
    qnull=function(x, p=0) -99
    doMethod = doMethod[doMethod!="Wassp1"]
  }  
  # Find out how many (non-chi-square) methods are implemented and what they are called

  TS_data=TS(x, (1:length(x))/(length(x)+1), phat(x), qnull)
  nummethods = length(TS_data)
  doMethodTS = names(TS_data)

  doMethodchi = c("EP large Pearson", "ES large Pearson", "EP small Pearson", "ES small Pearson",
        paste0(c("EP large LR", "ES large LR", "EP small LR", "ES small LR"),"-", ifelse(rate==0, "m", "p")) 
  )
  allMethods = c(doMethodTS, doMethodchi)  
  out = matrix(0, 2, nummethods+8) 
  colnames(out) = allMethods 

  if(doMethod[1]=="Default" || doMethod[1]=="all") {
    if(doMethod[1]=="Default") doMethod = c("W", "ZK", "ZC", "Wassp1", "EP small Pearson","ES small Pearson")      
    if(doMethod[1]=="all") doMethod = allMethods     
  } 

  if(any(doMethod %in% allMethods[1:nummethods]))
      out[ ,1:nummethods] = gof_cont(x, pnull, rnull, qnull, phat, TS, B)
  if(any(doMethod %in% allMethods[nummethods+1:8])) {  
      if(is.infinite(Range[1])) Range[1]=-99999
      if(is.infinite(Range[2])) Range[2]=99999
      out[ , nummethods+1:4] = t(chi_test_cont(x, pnull, phat(x),  formula="Pearson", 
                  nbins=nbins, Range=Range, Minimize=Minimize)[,c(1, 3)])
      out[ , nummethods+5:8] = t(chi_test_cont(x, pnull, phat(x),  formula="LR", rate=rate,
                    nbins=nbins, Range=Range, Minimize=Minimize)[,c(1, 3)])
  }  

  out = round(out[, doMethod, drop=FALSE], 4)
  list(Statistics=out[1, ], p.value=out[2, ])
}

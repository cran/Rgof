#' This function performs a number of gof tests for discrete data.
#' @param  x data set (the counts)
#' @param  pnull  cumulative distribution function under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  vals a vector of values of discrete random variables 
#' @param  phat  a function to estimate parameters from the data
#' @param  TS user supplied function to find test statistics
#' @param  nbins =c(100, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  B   =5000  number of simulation runs
#' @param  minexpcount =2 minimal expected bin count required
#' @param  maxProcessors =1 number of processors to use in parallel processing. If missing single processor is used.
#' @param  doMethod Methods to include in tests
#' @return A numeric matrix of test statistics and p values
#' @export
#' @examples 
#' # Tests to see whether data comes from a binomial (10, 0.5) distribution.
#' vals=0:10
#' pnull = function() pbinom(0:10, 10, 0.5)
#' rnull = function() table(c(0:10, rbinom(1000, 10, 0.5)))-1
#' x = rnull() 
#' gof_test_disc(x, pnull, rnull, vals, doMethod="all")
#' # Tests to see whether data comes from a binomial distribution with the success probability 
#' # estimated from the data.
#' pnull = function(p=0.5) pbinom(0:10, 10, ifelse(p>0&&p<1,p,0.001))
#' rnull = function(p=0.5) table(c(0:10, rbinom(1000, 10, ifelse(p>0&&p<1,p,0.001))))-1
#' phat = function(x) mean(0:10*x)/1000 
#' gof_test_disc(x, pnull, rnull, vals, phat)  

gof_test_disc <- function(x, pnull, rnull, vals, phat, TS,  
     nbins=c(100, 10), rate=0, B=5000, minexpcount=2.0, 
     maxProcessors=1, doMethod="Default") {

  NewTS = FALSE
  if(missing(TS)) TS = TS_disc
  else {
    NewTS=TRUE
    if(is.null(names(TS(x, (1:length(x))/length(x), nm_calc(sum(x)), vals)))) {
      message("result of TS has to be a named vector")
      return(NULL)
    } 
  }
  TS_data = TS(x, (1:length(x))/length(x), nm_calc(sum(x)), vals)
  nummethods = length(TS_data)
  doMethodTS = names(TS_data)
  doMethodchi = c("chi large Pearson", "chi small Pearson", 
      paste0(c("chi large LR", "chi small LR"),"-", ifelse(rate==0, "m", "p"))) 
  allMethods = c(doMethodTS, doMethodchi)
  out = matrix(0, 2, nummethods+4)    
  colnames(out) = allMethods
  if(doMethod[1]=="Default" || doMethod[1]=="all") {
    if(doMethod[1]=="Default") doMethod =c("K", "AD", "ZA", "ZC")
    if(doMethod[1]=="all") doMethod = allMethods      
  }  
  if(NewTS) {
    allMethods=doMethodTS
    doMethod=doMethodTS
  }
  
  if(!missing(phat)) check.functions(pnull, rnull, vals=vals, phat=phat, x=x)
  if(missing(phat)) check.functions(pnull, rnull, vals=vals, x=x)  
  
  Minimize = 1
  if(missing(phat)) {
    Minimize = 0
    phat = function(x) 0
  }

  if(any(doMethod %in% allMethods[1:nummethods])) {
     if(maxProcessors==1)
        out[ ,1:nummethods] = gof_disc(x, pnull, rnull, vals, phat, TS, rate, B)
     else {
        m = maxProcessors
        cl = parallel::makeCluster(m)
        z=parallel::clusterCall(cl, 
                            gof_disc, 
                            x = x,
                            pnull = pnull,
                            rnull = rnull,
                            vals = vals,  
                            phat = phat,
                            TS = TS,
                            rate=rate,
                            B = B/m
       )
       parallel::stopCluster(cl)
       #  Average power of cores
       tmp=0*z[[1]]
       for(i in 1:m) tmp=tmp+z[[i]]
       out[ ,1:nummethods] = tmp/m  
     }
  }

  if(any(doMethod %in% allMethods[nummethods+1:4])) {  
      out[ , nummethods+1:2] = t(chi_test_disc(x=x, pnull=pnull, phat(x), 
            nbins, "Pearson", rate, Minimize, minexpcount)[,c(1, 3)])
      out[ , nummethods+3:4] = t(chi_test_disc(x=x, pnull=pnull, phat(x), 
          nbins, "LR", rate, Minimize, minexpcount)[,c(1, 3)])      
  } 
  out = round(out[, doMethod, drop=FALSE], 4)
  list(Statistics=out[1, ], p.value=out[2, ])
}

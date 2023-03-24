#' This function performs a number of gof tests for discrete data.
#' @param  x data set (the counts)
#' @param  pnull  cumulative distribution function under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  vals a vector of values of discrete random variables 
#' @param  phat  a function to estimate parameters from the data
#' @param  nbins =c(100, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  B   =5000  number of simulation runs
#' @param  doMethod Methods to include in tests
#' @return A numeric matrix of test statistics and p values
#' @export
#' @examples 
#' # Tests to see whether data comes from a binomial (10, 0.5) distribution.
#' vals=0:10
#' pnull = function() pbinom(0:10, 10, 0.5)
#' rnull = function() table(c(0:10, rbinom(1000, 10, 0.5)))-1
#' x = rnull() 
#' gof_test_disc(x, pnull, rnull, vals)
#' # Tests to see whether data comes from a binomial distribution with the success probability 
#' # estimated from the data.
#' pnull = function(p=0.5) pbinom(0:10, 10, ifelse(p>0&&p<1,p,0.001))
#' rnull = function(p=0.5) table(c(0:10, rbinom(1000, 10, ifelse(p>0&&p<1,p,0.001))))-1
#' phat = function(x) mean(0:10*x)/1000 
#' gof_test_disc(x, pnull, rnull, vals, phat)  

gof_test_disc <- function(x, pnull, rnull, vals, phat,  
     nbins=c(100, 10), rate=0, B=5000, doMethod="Default") {

  doMethodTS = c("KS", "K","AD", "CvM",  "W", "ZA", "ZK",  "ZC", "Wassp1")    
  doMethodchi = c("chi large Pearson", "chi small Pearson", 
      paste0(c("chi large LR", "chi small LR"),"-", ifelse(rate==0, "m", "p"))) 
  allMethods = c(doMethodTS, doMethodchi)
  out = matrix(0, 2, 13)    
  colnames(out) = allMethods
  if(doMethod[1]=="Default" || doMethod[1]=="all") {
    if(doMethod[1]=="Default") doMethod =c("K", "AD", "ZA", "ZC")
    if(doMethod[1]=="all") doMethod = allMethods      
  }  

  if(!missing(phat)) check.functions(pnull, rnull, vals=vals, phat=phat, x=x)
  if(missing(phat)) check.functions(pnull, rnull, vals=vals, x=x)  
  
  Minimize = 1
  if(missing(phat)) {
    Minimize = 0
    phat = function(x) 0
  }

  if(any(doMethod %in% allMethods[1:9]))
      out[ ,1:9] = gof_disc(x=x, pnull=pnull, vals=vals, 
          rnull=rnull, phat=phat, rate=rate, B=B, doMethod=doMethodTS)
  if(any(doMethod %in% allMethods[10:13])) {  
      out[ , 10:11] = t(chi_test_disc(x=x, pnull=pnull, phat(x), 
            nbins, "Pearson", rate, Minimize)[,c(1, 3)])
      out[ , 12:13] = t(chi_test_disc(x=x, pnull=pnull, phat(x), 
          nbins, "LR", rate, Minimize)[,c(1, 3)])      
  } 
  doMethod = doMethod[doMethod!="ZK"]  
  out = round(out[, doMethod, drop=FALSE], 4)
  list(Statistics=out[1, ], p.value=out[2, ])
}

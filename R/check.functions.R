#' This function checks whether the inputs have the correct format
#' 
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  qnull  routine to calculate quantiles under null hypothesis
#' @param  phat   function to estimate parameters from the data
#' @param  vals   vector of discrete values
#' @param  x      data
#' @return NULL
#' 
check.functions=function(pnull, rnull, qnull, phat, vals, x) {
  npnull = length(formals(pnull))
  nrnull = length(formals(rnull))
# Continuous -  No Estimation 
  if(missing(phat) && missing(vals)) {
    if(npnull!=1) message("pnull should have one argument for simple hypothesis if data is continuous")
    if(nrnull!=0) message("rnull should have no argument for simple hypothesis if data is continuous")
    if(!missing(qnull)) {
      nqnull = length(formals(qnull))
      if(nqnull!=1) message("qnull should have one argument for simple hypothesis")
    }  
  } 
# Continuous -  With Estimation 
  if(!missing(phat) && missing(vals)) {
    if(npnull!=2) message("pnull should have two arguments for composite hypothesis if data is continuous")
    if(nrnull!=1) message("rnull should have one argument for composite hypothesis if data is continuous")
    if(!missing(qnull)) {
      nqnull = length(formals(qnull))
      if(nqnull!=2) message("qnull should have two arguments for composite hypothesis")
    }
    if(is.function(phat) && length(formals(phat))!=1) 
      message("phat should have one argument x, the data")
  }  
  
# Discrete -  No Estimation 
  if(missing(phat) && !missing(vals)) {
    if(length(vals)!=length(pnull())) message("vals and pnull() should have the same length")
    if(length(vals)!=length(rnull())) message("vals and rnull() should have the same length")
    if(npnull!=0) message("pnull should have no argument for simple hypothesis if data is discrete")
    if(nrnull!=0) message("rnull should have no argument for simple hypothesis if data is discrete")
    if(!missing(x) && length(x)!=length(vals))
      message("For discrete data x and vals have to have the same length") 
    if(min(diff(pnull()))<0) message("For discrete data pnull should be strictly increasing")
    if(min(diff(vals))<0) message("For discrete data vals should be strictly increasing")
  } 
# Discrete -  With Estimation 
  if(!missing(phat) && !missing(vals)) {
    if(npnull!=1) message("pnull should have one arguments for composite hypothesis if data is discrete")
    if(nrnull!=1) message("rnull should have one argument for composite hypothesis if data is discrete")
    if(is.function(phat) && length(formals(phat))!=1) 
        message("phat should have one argument x, the data")
    if(!missing(x) && length(x)!=length(vals))
      message("For discrete data x and vals have to have the same length") 
    if(min(diff(vals))<0) message("For discrete data vals should be strictly increasing")
  }  
}

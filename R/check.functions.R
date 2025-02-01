#' This function checks whether the inputs have the correct format
#' 
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  phat   =function(x) -99, function to estimate parameters from the data, or -99
#' @param  vals   vector of discrete values
#' @param  x      data
#' @return NULL
#' 
check.functions=function(pnull, rnull, phat=function(x) -99, vals, x) {
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if(all(is.wholenumber(x)) & missing(vals)) 
     message("vector x is all integers but argument vals is missing. Is data discrete?")
  if(!all(is.wholenumber(x)) & !missing(vals)) 
    message("vector x has to be all integers/counts for discrete data")
  npnull = length(formals(pnull))
  nrnull = length(formals(rnull))
# Continuous -  No Estimation 
  if(phat(x)[1]==-99 && missing(vals)) {
    if(npnull!=1) message("pnull should have one argument for simple hypothesis if data is continuous")
    if(nrnull!=0) message("rnull should have no argument for simple hypothesis if data is continuous")
  } 
# Continuous -  With Estimation 
  if(phat(x)[1]!=-99 && missing(vals)) {
    if(npnull!=2) message("pnull should have two arguments for composite hypothesis if data is continuous")
    if(nrnull!=1) message("rnull should have one argument for composite hypothesis if data is continuous")
    if(is.function(phat) && length(formals(phat))!=1) 
      message("phat should have one argument x, the data")
  }  
# Discrete -  No Estimation 
  if(phat(x)[1]==-99 && !missing(vals)) {
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
  if(phat(x)[1]!=-99 && !missing(vals)) {
    if(npnull!=1) message("pnull should have one argument for composite hypothesis if data is discrete")
    if(nrnull!=1) message("rnull should have one argument for composite hypothesis if data is discrete")
    if(is.function(phat) && length(formals(phat))!=1) 
        message("phat should have one argument x, the data")
    if(!missing(x) && length(x)!=length(vals))
      message("For discrete data x and vals have to have the same length") 
    if(min(diff(vals))<0) message("For discrete data vals should be strictly increasing")
  }  
  
}

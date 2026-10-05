#' MLE of ratio of mixture distributions
#' 
#' This function uses the Newton-Raphson algorithm
#' to find the mle of a mixture distribution
#'
#' @param fx density 1 evaluated at data set
#' @param gx density 2 evaluated at data set
#' @keywords internal
#' @return a number
#' @export 
mlemix=function(fx, gx, acc=0.001) {
  ow=0.5
  k=0
  repeat {
    k=k+1
    tmp=(fx-gx)/(ow*fx+(1-ow)*gx)
    nw=ow+sum(tmp)/sum(tmp^2)
    if(abs(ow-nw)< acc) break
    if(k>100) break
    ow=nw
  }
  pmin(1, pmax(0, nw))
}
#' mle's of truncated exponential distribution
#' 
#' This function uses the Newton-Rahson algorithm
#' to find the mle's of a truncated exponential distribution
#'
#' @param x data set
#' @keywords internal
#' @return a vector
#' @export 
mletexp=function(x) {
  mx=mean(x)
  Mx=max(x)*(1+1/length(x))
  ow=1/mx 
  k=0
  repeat {
    k=k+1
    z=exp(-ow*Mx)
    den=1/ow-mx-Mx*z/(1-z)
    num=-1/ow^2+Mx^2*z/(1-z)^2
    nw=ow-den/num
    if(abs(ow-nw)<0.001 | k>100) break
    ow=nw
  }
  c(ifelse(nw<0.01, 0.01, nw), Mx)
}
#' mle of double truncated exponential distribution
#' 
#' This function uses the Newton-Rahson algorithm
#' to find the mle of a doubly truncated exponential distribution
#'
#' @param x data set
#' @keywords internal
#' @return a vector
#' @export 
mledexp=function(x) {
  mx=mean(x)
  ow=1/mx 
  k=0
  repeat {
    k=k+1
    u=exp(-ow)
    v=exp(-2*ow)
    den=1/ow-mx+(u-2*v)/(u-v)
    num=-1/ow^2+u*v/(u-v)^2
    nw=ow-den/num
    if(abs(ow-nw)<0.001 | k>100) break
    ow=nw
  }
  ifelse(nw<0.01, 0.01, nw)
}

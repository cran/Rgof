#' Testing routine for continuous data
#' 
#' find test statistic or p value for continuous data
#'
#' @param  x data set 
#' @param  pnull function to find cdf under  null hypothesis
#' @param p =0, parameter estimates or 0
#' @param TSextra list with additional info
#' @return a number, value fo statistic or p value
#' @export
myTS_cont=function(x, pnull, p=0, TSextra) {
  nbins=TSextra$nbins
  out=matrix(0, 2, length(nbins))
  colnames(out)=paste0("Chi_C", nbins)
  for(i in seq_along(nbins)) {
    bins=quantile(x, c(0:nbins[i])/nbins[i])
    bins[1]=-Inf;bins[nbins[i]+1]=Inf
    O=hist(x, bins, plot=FALSE)$counts
    if(length(formals(pnull))==1)
      E=length(x)*diff(pnull(bins))
    else E=length(x)*diff(pnull(bins, p))
    out[1, i]=sum((O-E)^2/E)
    if(length(formals(pnull))==1) m=0
    else m=length(p)
    out[2, i]=1-pchisq(out[1, i], length(O)-1-m)
  }  
  if(TSextra$statistic) return(round(out[1, ], 4))
  round(out[2, ],4)
}
#' Testing routine for discrete data
#' 
#' find test statistic or p value for discrete data
#'
#' @param  x data set 
#' @param  pnull function to find cdf under  null hypothesis
#' @param  p =0, parameter estimates or 0
#' @param  TSextra list with additional info
#' @param  vals =NA values of discrete random variable, or NA
#' @return a number, value fo statistic or p value
#' @export
myTS_disc=function(x, pnull, p=0, vals, TSextra) {
  O=x
  if(length(formals(pnull))==0)
    E=sum(x)*diff(c(0, pnull()))
  else E=sum(x)*diff(c(0, pnull(p)))
  chi=sum((O-E)^2/E)
  names(chi)="Chi_D"
  if(TSextra$statistic) return(round(chi, 4))
  if(length(formals(pnull))==0) m=0
  else m=length(p)
  pval=1-pchisq(chi, length(O)-1-m)
  names(pval)="Chi_D"
  round(pval, 4)
}

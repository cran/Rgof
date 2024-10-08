#' This function creates several type of bins for continuous data
#' @param  x data set
#' @param  pnull  cdf under the null hypothesis
#' @param  qnull =NA quantile function, if available
#' @param  phat =function(x) -99  parameters for pnull
#' @param  DataBased =FALSE bins based on data, not expected counts
#' @param  nbins =c(50, 10) number of bins
#' @param  minexpcount =5 smallest expected count per bin
#' @param  Range  =c(-99999, 99999) limits of possible observations, if any
#' @return A list of bins and bin probabilities

make_bins_cont = function(x, pnull, qnull=NA, phat=function(x) -99, 
                          DataBased=FALSE, nbins=c(50,10), 
                          minexpcount=5, Range=c(-99999, 99999)) {
                          
  n = length(x)
  res = formals(pnull)
  param = phat(x)
  if(length(res)==1) pf=function(x) pnull(x)
  else pf=function(x) pnull(x, param)
  if(is.function(qnull)) {
    if(length(res)==1) qf=function(x) qnull(x)
    else qf=function(x) qnull(x, param)
  }
  out = as.list(1:4)
  names(out) = c("ES-l", "ES-s", "EP-l", "EP-s")
  combine.bins = function(E, bins, minexpcount=5, maxbins=1000) {
    n = sum(E)
    k = length(E)
    while ( (min(E)<minexpcount) || (k>maxbins) ){
      i = which.min(E)[1]
      if(i==1) {
        E = c(E[1]+E[2], E[3:k])
        bins = bins[-2]
      }
      if(i==2) {
        if((E)[1]<(E)[3]) {
          E = c(E[1]+E[2], E[3:k])
          bins=bins[-2]
        }
        else {
          E = c(E[1], E[2]+E[3], E[4:k])
          bins=bins[-3]
        }
      }
      if(i==k-1) {
        if(E[k-2]<E[k]) {
          E = c(E[1:(k-3)], E[k-2]+E[k-1], E[k])
          bins=bins[-(k-1)]
        }
        else {
          E = c(E[1:(k-2)], E[k-1]+E[k])
          bins=bins[-k]
        }
      }
      if(i==k) {
        E = c(E[1:(k-2)], E[k-1]+E[k])
        bins=bins[-k]
        
      }
      if(i>2 && i<k-1) {
        if( (E)[i-1]<(E)[i+1] ) {
          E = c(E[1:(i-2)], E[i-1]+E[i], E[(i+1):k])
          bins=bins[-i]
        }
        else {
          E = c(E[1:(i-1)], E[i]+E[i+1], E[(i+2):k])
          bins=bins[-(i+1)]
        }
      }
      k = k-1
    }
    bins
  }

  for(i in 1:2) {
    if(Range[1]==-99999&Range[2]==99999) {
      bins = seq(min(x)-1e-05, max(x)+1e-05, length=nbins[i]+1)
      bins[c(1, nbins[i]+1)] = c(-Inf, Inf)
    }  
    if(Range[1]!=-99999&Range[2]==99999) {
      bins = seq(Range[1], max(x)+1e-05, length=nbins[i]+1)   
      bins[nbins[i]+1] = Inf
    }  
    if(Range[1]==-99999&Range[2]!=99999) {
      bins = seq(min(x)-1e-05, Range[2], length=nbins[i]+1)
      bins[1] = Inf 
    }  
    if(Range[1]!=-99999&Range[2]!=99999)
      bins = seq(Range[1], Range[2], length=nbins[i]+1)
    if(!DataBased) E=n*diff(pf(bins))
    else E=bincounter(x, bins)
    out[[i]] = combine.bins(E, bins, minexpcount, nbins[i])
  } 
  
  for(i in 1:2) {
    if(!is.function(qnull)) bins = c(-Inf, quantile(x, 1:(nbins[i]-1)/nbins[i]), Inf)
    else bins = qf(0:nbins[i]/nbins[i])
    if(Range[1]!=-99999) bins[1]=Range[1]
    if(Range[2]!=99999) bins[nbins[i]+1]=Range[2]
    if(!DataBased) E=n*diff(pf(bins))
    else E=bincounter(x, bins)
    out[[i+2]] = combine.bins(E, bins, minexpcount, nbins[i])
  }
  out
}

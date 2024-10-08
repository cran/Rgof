#' This function creates several types of bins for discrete data
#' @param  x counts
#' @param  pnull cumulative distribution function 
#' @param  phat =function(x) -99, function to estimated parameters, or -99
#' @param  nbins =c(50, 10) number of bins
#' @param  minexpcount =5 smallest expected count per bin
#' @return A list of indices

make_bins_disc = function(x, pnull, phat=function(x) -99, nbins=c(50,10), minexpcount=5) {
  n = sum(x)
  out = as.list(1:2)
  names(out) = c("l", "s")
  res = formals(pnull)
  if(length(res)==0) pf=pnull()
  else pf=pnull(phat(x))
  combine.bins = function(E, minexpcount=5, maxbins=1000) {
    k = length(E)
    bins = as.list(1:k)
    while (min(E)<minexpcount){
      i = which.min(E)[1]
      if(length(E)!=k) break
      if(any(is.na(E))) return(0)
      if(i==1) {
        E = c(E[1]+E[2], E[3:k])
        bins[[1]] = c(bins[[1]], bins[[2]])
        bins = bins[-2]
      }
      if(i==k) {
        E = c(E[1:(k-2)], E[k-1]+E[k])
        bins[[k]] = c(bins[[k-1]], bins[[k]])
        bins = bins[-(k-1)]
      }
      if(i>1 && i<k) {
        if( (E)[i-1]<(E)[i+1] ) {
          if(i==2) E = c(E[1]+E[2], E[3:k])
          else E = c(E[1:(i-2)], E[i-1]+E[i], E[(i+1):k])
          bins[[i]] = c(bins[[i-1]], bins[[i]])
          bins = bins[-(i-1)]
        }
        else {
          if(i==k-1) E = c(E[1:(k-2)], E[k-1]+E[k])
          else E = c(E[1:(i-1)], E[i]+E[i+1], E[(i+2):k])
          bins[[i]] = c(bins[[i]], bins[[i+1]])
          bins = bins[-(i+1)]
        }
      }
      k = k-1
    }
    while (k>maxbins) {
       binlengths=unlist(lapply(bins, length))
       shortestbins=min(binlengths)
       Index=c(1:k)[binlengths==min(binlengths)]
       if(length(Index)==1) i=Index
       else i=sample(Index,1)
       if(length(E)!=k) break
       if(any(is.na(E))) return(0)
       if(i==1) {
         E = c(E[1]+E[2], E[3:k])
         bins[[1]] = c(bins[[1]], bins[[2]])
         bins = bins[-2]
       }
       if(i==k) {
         E = c(E[1:(k-2)], E[k-1]+E[k])
         bins[[k]] = c(bins[[k-1]], bins[[k]])
         bins = bins[-(k-1)]
       }
       if(i>1 && i<k) {
         if( (E)[i-1]<(E)[i+1] ) {
           if(i==2) E = c(E[1]+E[2], E[3:k])
           else E = c(E[1:(i-2)], E[i-1]+E[i], E[(i+1):k])
           bins[[i]] = c(bins[[i-1]], bins[[i]])
           bins = bins[-(i-1)]
         }
         else {
           if(i==k-1) E = c(E[1:(k-2)], E[k-1]+E[k])
           else E = c(E[1:(i-1)], E[i]+E[i+1], E[(i+2):k])
           bins[[i]] = c(bins[[i]], bins[[i+1]])
           bins = bins[-(i+1)]
         }
       }
       k = k-1
    }
    bins
  }
  for(i in 1:2) {
    out[[i]] = combine.bins(n*diff(c(0,pf)), minexpcount, nbins[i])
  } 
  out
}

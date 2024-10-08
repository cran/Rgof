#' This function performs a number of chi-square gof tests for continuous data
#' @param  x data set
#' @param  pnull  cdf under the null hypothesis
#' @param  w  function to find weights of observations, returns -99 if data is unweighted
#' @param  phat =function(x) -99, estimated parameters, or starting values of multi-D minimum chi square minimization, or -99 if no estimation is done
#' @param  qnull =NA quantile function, if available
#' @param  nbins =c(50, 10) number of bins for chi-square tests
#' @param  rate =0, rate of Poisson if sample size is random
#' @param  Range  =c(-99999, 99999) limits of possible observations, if any
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat =TRUE, if TRUE param is estimated parameters and no minimization is used
#' @param  allbins set of bins to use
#' @return A numeric matrix of test statistics, degrees of freedom and p.values

chi_test_cont = function(x, pnull, w=function(x) -99, phat=function(x) -99, 
                         qnull=NA, nbins=c(50,10), rate=0, 
                         Range=c(-99999, 99999), minexpcount=5, 
                         ChiUsePhat=TRUE, allbins) {
  out = matrix(0, 8, 3)
  colnames(out) = c("Statistic",  "p.values", "df")
  rownames(out) = c("ES-l-P", "ES-s-P", "EP-l-P", "EP-s-P",
                    "ES-l-L", "ES-s-L", "EP-l-L", "EP-s-L")
  if(length(formals(w))==1) wx=w(x)
  else wx=w(x, phat(x)) # with parameter estimation  
  if(missing(allbins)) {
     allbins = make_bins_cont(x, pnull, 
              qnull, phat, 
              DataBased=ifelse(wx[1]==-99, FALSE, TRUE), 
              nbins, minexpcount, Range)
  } 
  chi_stat = function(param, O, bins, pnull, w, formula="P", rate=0) {
    n = sum(O)
    k = length(bins)-1
    res = formals(pnull)
    if(length(res)==1) p=pnull(bins)
    else p=pnull(bins, param)
    if(rate==0) rate=n
    if(length(formals(w))==1) wx=w(x)
    else wx=w(x, param) 
    if(wx[1]==-99) E = rate*diff(p)
    else {
        xrange = diff(range(bins[is.finite(bins)]))
        h = 0.001*xrange
        E=rep(0, k)
        if(length(res)==1) f=function(x) (pnull(x+h)-pnull(x))/h/w(x)
        else f=function(x) (pnull(x+h, param)-pnull(x, param))/h/w(x,param)
        if(is.finite(bins[1])) 
           E[1]=rate*stats::integrate(f, bins[1], bins[2])$value
        else 
           E[1]=rate*stats::integrate(f, min(x)-xrange/10, bins[2])$value
        if(is.finite(bins[k+1])) 
          E[k]=rate*stats::integrate(f, bins[k], bins[k+1])$value
        else 
          E[k]=rate*stats::integrate(f, bins[k], max(x)+xrange/10)$value
        for(i in 2:(k-1)) E[i]=rate*stats::integrate(f, bins[i], bins[i+1])$value
    }
    if(formula=="P") chi = sum((O-E)^2/E)
    else {
      I = seq_along(O)[O>0]
      chi = 2*sum(E[I]-O[I]+O[I]*log(O[I]/E[I]))
    }
    chi
  }
  formula=c("P", "L")
  type=c("ES", "EP")
  binsize=c("l", "s")  
  for(i in 1:2) {
    for(j in 1:2) {
      for(k in 1:2) {
        bins  = allbins[[paste0(type[j],"-",binsize[k])]]
        N=bincounter(x, bins)
        if(ChiUsePhat) 
            chi=chi_stat(phat(x), N, bins, pnull, w, formula[i], rate)
        else chi=stats::optim(phat(x), chi_stat, method = "BFGS", 
                       O=N, bins=bins, pnull=pnull, w=w,
                       formula=formula[i], rate=rate)$value    
        adj.df=1     
        df = length(N) -  
             ifelse(rate==0, 1, 0) - 
             ifelse(phat(x)[1]==-99, 0, length(phat(x))) 
        out[paste0(type[j],"-",binsize[k],"-",formula[i]), ] =  c(chi, 1-stats::pchisq(chi, df), df)
        }
    }
  }

  out
}

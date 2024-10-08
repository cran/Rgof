#' This function performs a number of chi-square gof tests for continuous data
#' @param  x data set
#' @param  pnull  cdf under the null hypothesis
#' @param  phat  =function(x) -99, function to estimate parameters, or starting values of multi-D minimum chi square minimization, or -99 if no parameters are estimated
#' @param  nbins =c(50, 10) number of bins for chi-square tests
#' @param  rate =0, rate of Poisson if sample size is random
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  allbins set of bins to use
#' @return A numeric matrix of test statistics, degrees of freedom and p.values

chi_test_disc = function(x, pnull, phat=function(x) -99, nbins=c(50,10), 
        rate=0, minexpcount=5, ChiUsePhat=TRUE, allbins) {
  out = matrix(0, 4, 3)
  colnames(out) = c("Statistic",  "p.values", "df")
  rownames(out) = c("l-P", "s-P", "l-L", "s-L")
  if(missing(allbins)) allbins = make_bins_disc(x, pnull, phat, nbins, minexpcount)
  chi_stat = function(param, x, pnull, bins, formula="P", rate=0) {
      n = sum(x)
      if(rate==0) rate=n
      res = formals(pnull)
      if(length(res)==0) p=pnull()
      else p=pnull(param)
      O  = unlist(lapply(bins, function(q) sum(x[q])))
      E = n*diff(c(0,p))
      E = unlist(lapply(bins, function(q) sum(E[q])))
      if(formula=="P") chi = sum((O-E)^2/E)
      else {
        I = seq_along(O)[O>0]
        chi = 2*sum(E[I]-O[I]+O[I]*log(O[I]/E[I]))
      } 
      chi
  }
  formula=c("P", "L")
  binsize=c("l", "s")  
  for(i in 1:2) {
    for(k in 1:2) {
        bins  = allbins[[binsize[k]]]
        if(ChiUsePhat) 
           chi=chi_stat(phat(x), x, pnull, bins, formula[i], rate)
        else chi=stats::optim(phat(x), chi_stat, method = "BFGS", 
             x=x, pnull=pnull, bins=bins,
             formula=formula[i], rate=rate)$value
        df = length(bins) - 
            ifelse(phat(x)[1]==-99, 0, length(phat(x))) - 
            ifelse(rate==0, 1, 0)
        out[paste0(binsize[k],"-",formula[i]), ] =  c(chi, 1-stats::pchisq(chi, df), df)
    }
  }
  out
}

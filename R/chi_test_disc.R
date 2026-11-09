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
#' @noRd

chi_test_disc = function(x, pnull, phat=function(x) -99, nbins=c(50,10), 
        rate=0, minexpcount=5, ChiUsePhat=TRUE, allbins) {
  out = matrix(0, 4, 3)
  param0 = phat(x)
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
        if(ChiUsePhat || param0[1]==-99) {
           chi=chi_stat(param0, x, pnull, bins, formula[i], rate)
        } else {
           f = function(param) {
             z = tryCatch(chi_stat(param, x, pnull, bins, formula[i], rate),
                          error=function(e) Inf, warning=function(e) Inf)
             if(length(z)!=1 || !is.finite(z) || z<0) Inf else z
           }
           chi0=f(param0)
           chi=chi0
           if(length(param0)==1) {
             h=0.1*max(1, abs(param0))
             a=param0-h; b=param0+h
             for(ii in 1:12) {
               fa=f(a); fb=f(b)
               if(is.finite(fa) && is.finite(fb) && fa>=chi0 && fb>=chi0) break
               h=2*h; a=param0-h; b=param0+h
             }
             # f is Inf where the model can not be evaluated, optimize warns about that
             fit=tryCatch(suppressWarnings(stats::optimize(f, interval=c(a,b))), error=function(e) NULL)
             if(!is.null(fit) && is.finite(fit$objective) && fit$objective<=chi0) chi=fit$objective
           } else {
             fit=tryCatch(stats::optim(param0, f, method="Nelder-Mead"), error=function(e) NULL)
             if(!is.null(fit) && is.finite(fit$value) && fit$value<=chi0) chi=fit$value
           }
        }
        df = length(bins) - 
            ifelse(param0[1]==-99, 0, length(param0)) - 
            ifelse(rate==0, 1, 0)
        out[paste0(binsize[k],"-",formula[i]), ] =  c(chi, 1-stats::pchisq(chi, df), df)
    }
  }
  out
}

#' This function estimates the power of test routines that calculate p value(s)
#' @param TS routine to calculate test statistics.
#' @param vals =NA if data is discrete, a vector of possible values  
#' @param pnull routine to calculate the cdf under the null hypothesis
#' @param ralt  generate data under altenative hypothesis
#' @param param_alt values of parameter under the alternative hypothesis. 
#' @param phat function to estimate parameters, function(x) -99 if no parameter estimation
#' @param TSextra list (possibly) passed to TS
#' @param alpha =0.05  type I error.
#' @param B = 1000 number of simulation runs to estimate the power.
#' @return A matrix of power values

power_newtest <- function(TS, vals=NA, pnull, ralt, param_alt, phat, TSextra, alpha=0.05, B=1000) {
     Continuous=ifelse(is.na(vals)[1], TRUE, FALSE)
     WithTSextra=ifelse("aaa"%in%names(TSextra), FALSE, TRUE)
     dta=ralt(param_alt[1])
     if(Continuous && !WithTSextra) tmp=TS(dta, pnull, phat(dta)) 
     if(Continuous && WithTSextra) tmp=TS(dta, pnull, phat(dta), TSextra) 
     if(!Continuous && !WithTSextra) tmp=TS(dta, pnull, phat(dta), vals) 
     if(!Continuous && WithTSextra) tmp=TS(dta, pnull, phat(dta), vals, TSextra) 
     pwr=matrix(0, length(param_alt), length(TS))
     rownames(pwr)=param_alt
     colnames(pwr)=names(tmp)
     for(j in 1:length(param_alt)) {
        for(i in 1:B) {
          dta=ralt(param_alt[j])
          if(Continuous && !WithTSextra) tmp=TS(dta, pnull, phat(dta)) 
          if(Continuous && WithTSextra) tmp=TS(dta, pnull, phat(dta), TSextra) 
          if(!Continuous && !WithTSextra) tmp=TS(dta, pnull, phat(dta), vals) 
          if(!Continuous && WithTSextra) tmp=TS(dta, pnull, phat(dta), vals, TSextra) 
          pwr[j, ]=pwr[j, ]+ifelse(tmp<alpha,1,0) 
        }
     }
     pwr/B
}
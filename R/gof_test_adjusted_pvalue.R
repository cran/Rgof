#' Adjusted p values for simultaneous testing in the goodness-of-fit problem.
#' 
#' This function performs a number of goodness-of-fit tests and finds the adjusted p value for the combined test.
#' 
#' For details on the usage of this routine consult the vignette with vignette("Rgof","Rgof")
#' 
#' @param  x data set
#' @param  vals =NA, values of discrete RV, or NA if data is continuous
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  w (Optional) function to calculate weights, returns -99 if no weights
#' @param  phat  =function(x) -99, function to estimate parameters from the data, or -99 if no parameters are estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  TSextra =NA, list passed to TS, if desired, or NA
#' @param  nbins =c(100, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any, for chi-square tests
#' @param  B   =c(5000,1000)  number of simulation runs for individual and for adjusted p values
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  maxProcessor number of cores to use
#' @param  doMethods a vector of codes for the methods to include. If missing, a default selection of methods are used.
#' @return None 
#' @export
#' @examples
#' # Tests to see whether data comes from a standard normal distribution.
#' pnull = function(x) pnorm(x)
#' rnull = function()  rnorm(100)
#' x = rnorm(100)
#' gof_test_adjusted_pvalue(x, NA, pnull, rnull, B=c(500, 200), 
#'                             maxProcessor=1)
#' # Tests to see whether data comes from a normal distribution with standard deviation 1 
#' # and the mean estimated.
#' pnull=function(x, m) pnorm(x, m)
#' rnull=function(m) rnorm(100, m)
#' TSextra = list(qnull=function(x, m=0) qnorm(x, m))
#' phat=function(x) mean(x)
#' x = rnorm(100, 1, 2)
#' gof_test_adjusted_pvalue(x, NA, pnull, rnull, phat=phat, 
#'                         TSextra=TSextra, B=c(500, 200), maxProcessor=1)
#' # Tests to see whether data comes from a binomial (10, 0.5) distribution.
#' vals=0:10
#' pnull = function() pbinom(0:10, 10, 0.5)
#' rnull = function() table(c(0:10, rbinom(1000, 10, 0.5)))-1
#' x = rnull() 
#' gof_test_adjusted_pvalue(x, vals, pnull, rnull, 
#'                         B=c(500, 200), maxProcessor=1)
#' # Tests to see whether data comes from a binomial distribution with 
#' # the success probability estimated from the data.
#' pnull = function(p=0.5) pbinom(0:10, 10, p)
#' rnull = function(p=0.5) table(c(0:10, rbinom(1000, 10, p)))-1
#' phat=function(x) mean(rep(0:10,x))/10 
#' gof_test_adjusted_pvalue(x, vals, pnull, rnull, phat=phat, 
#'                         B=c(500, 200), maxProcessor=1) 
gof_test_adjusted_pvalue <- function(x, vals= NA, pnull, rnull, 
                    w=function(x) -99, phat=function(x) -99, 
                    TS, TSextra=NA, nbins=c(50, 10), rate=0, 
                    Range=c(-Inf, Inf), B=c(5000,1000),  minexpcount=5.0,  
                    ChiUsePhat=TRUE, maxProcessor, doMethods) {
 
   if(length(B)==1) B=c(B, B) # this routine needs two simulation sizes
   if(any(is.na(vals))) {
     Continuous=TRUE
     dta=list(x=x)
     check.functions(pnull, rnull, phat, x=x)
   }  
   else {
     Continuous=FALSE
     dta=list(x=x, vals=vals)
     check.functions(pnull, rnull, phat, vals, x)
   }  
   if(any(is.na(TSextra))) TSextra=list(pnull=pnull, phat=phat, 
                                    w=w, Continuous=Continuous)
   else TSextra = c(TSextra, pnull=pnull, phat=phat, 
                    w=w, Continuous=Continuous)
   Noqnull=FALSE
   if(!("qnull" %in% names(TSextra))) {
     Noqnull=TRUE
     TSextra=c(TSextra, qnull=function(x) -99)
   }
   Noqnull=FALSE
   if(!("qnull" %in% names(TSextra))) {
     Noqnull=TRUE
     TSextra=c(TSextra, qnull=function(x) -99)
   }   
   WithWeights = TRUE
   if(length(formals(w))==1) {
     if(w(x[1])==-99) WithWeights = FALSE
   }
   if(length(x)>10000 && maxProcessor==1)
     message("Consider using parallel processing with maxProcessor= (your number of cores)")
   # adjust number of bins to account for parameter estimation
   if(abs(phat(x)[1]+99)<0.001) nbins=nbins+length(phat(x))
   if(missing(TS)) {
     if(Continuous) {
       if(!WithWeights) { #data is not weighted
         typeTS=1
         TS = TS_cont
         if(missing(doMethods)) doMethods=c("W", "ZC", "AD", "ES-s-P")
       }
       else {
         typeTS=2
         TS = TSw_cont
       }
     }
     else {
       typeTS = 5
       TS = TS_disc
       if(missing(doMethods)) doMethods=c("W", "AD", "s-P")
     }
   }   
   else {
     # can't do parallel processing if TS written in C/C++
     if(substr(deparse(TS)[2], 1, 5)==".Call") {
       message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
       maxProcessor=1
     }
     if(Continuous) {
       if(length(formals(TS))>4) {
         message("TS for continuous data should have either 3 or 4 arguments")
         return(NULL)
       }  
       typeTS=length(formals(TS))
     }  
     else {
       if(length(formals(TS))>6) {
         message("TS for discrete datashould have either 4 or 5 arguments")
         return(NULL)
       }
       typeTS=length(formals(TS))+1
     }
   }
   TS_data=calcTS(dta, TS, typeTS, TSextra)  
   if(is.null(names(TS_data))) {
     message("result of TS has to be a named vector")
     return(NULL)
   }
   if(missing(maxProcessor)) {
     maxProcessor=parallel::detectCores(logical = FALSE)-1
     message(paste("Using",maxProcessor,"cores ..."))
   } 
   if(maxProcessor==1) {
      A=gof_adj_C1(dta, rnull, vals, TS, typeTS, TSextra, B[1])
      pvals=gof_adj_C2(dta, rnull, vals, TS, typeTS, TSextra, A, B[2])
   }    
   else {
     cl = parallel::makeCluster(maxProcessor)
     z=parallel::clusterCall(cl, gof_adj_C1, 
                dta, rnull, vals, TS, typeTS, TSextra, 
                             B=round(B[1]/maxProcessor)
     )
     A=z[[1]]
     for(i in 1:maxProcessor) A=rbind(A, z[[i]])
     z=parallel::clusterCall(cl, gof_adj_C2, 
                          dta, rnull, vals, TS, typeTS, TSextra, 
                          A, B=round(B[2]/maxProcessor)
     )
     parallel::stopCluster(cl)
     pvals=z[[1]]
     for(i in 1:maxProcessor) pvals=rbind(pvals,z[[i]])
     pvals=pvals
   }
   colnames(pvals)=names(TS_data)
   if(typeTS==1 | typeTS==5) { 
     B[2]=nrow(pvals)-1
     if(typeTS==1) chipvals=matrix(0, B[2]+1, 8)
     if(typeTS==5) chipvals=matrix(0, B[2]+1, 4)            
     for(i in 1:(B[2]+1)) {
        if(i==1) xsim=x
        else {
          if(length(formals(rnull))==0) xsim=rnull()
          else xsim=rnull(TSextra$phat(x))
          if(TSextra$Continuous) dtasim=list(x=xsim)
          else dtasim=list(x=xsim, vals=vals)
        }      
        if(typeTS==1) {
          if(is.infinite(Range[1])) Range[1]=-99999
          if(is.infinite(Range[2])) Range[2]=99999
          tmp = chi_test_cont(xsim, TSextra$pnull, TSextra$w, TSextra$phat,  
                           ifelse(TSextra$Noqnull, NA, TSextra$qnull),
                           nbins, rate, Range, minexpcount, ChiUsePhat)
          chipvals[i, ] = tmp[, 2]
        }
        if(typeTS==5) {
          tmp=chi_test_disc(xsim, TSextra$pnull, TSextra$phat, 
                         nbins, rate, minexpcount, ChiUsePhat)
          chipvals[i, ] = tmp[, 2]
        }                  
     }
      colnames(chipvals)=rownames(tmp)
      pvals=cbind(pvals, chipvals, 4)
   }
   if(missing(doMethods)) doMethods=colnames(pvals)
   if(test_methods(doMethods, Continuous, WithWeights)) return(NULL)
   pvals = round(pvals[, doMethods, drop=FALSE],4)
   minp_x=min(pvals[1, ])
   minp_sim=apply(pvals[-1, ,drop=FALSE], 1, min)
   z=seq(0, 1, length=250)
   y=z
   for(i in 1:250) y[i]=sum(minp_sim<=z[i])/B[2]
   I=c(1:250)[z>minp_x][1]-1
   slope=(y[I+1]-y[I])/(z[I+1]-z[I])
   minp_adj=round(y[I]+slope*(minp_x-z[I]),4)
   message("p values of individual tests:")
   for(i in 1:ncol(pvals)) message(paste(doMethods[i],": ", pvals[1,i]))
   message(paste0("adjusted p value of combined tests: ", minp_adj))
   
}

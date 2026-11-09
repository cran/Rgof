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
#' @param  TSextra list passed to TS, if desired, or missing
#' @param  qnull (optional) quantile function under the null hypothesis, with the same arguments as
#'   pnull. Needed for the Wasserstein test (Wassp1) with continuous data; if it is missing that test
#'   is not run. It can also be given as TSextra=list(qnull=...), as in earlier versions.
#' @param  nbins =c(100, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any, for chi-square tests
#' @param  B   =c(5000,1000)  number of simulation runs for individual and for adjusted p values
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = FALSE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  maxProcessor number of cores to use
#' @param  SuppressMessages =FALSE
#' @param  seed =NULL, optional seed for reproducibility
#' @param  doMethods a vector of codes for the methods to include. If missing, a default selection of methods are used.
#' @return (invisibly) a list with the p values of the individual tests and the adjusted p value. These are also printed as messages.
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
#' qnull = function(x, m=0) qnorm(x, m)
#' phat=function(x) mean(x)
#' x = rnorm(100, 1, 2)
#' gof_test_adjusted_pvalue(x, NA, pnull, rnull, phat=phat, 
#'                         qnull=qnull, B=c(500, 200), maxProcessor=1)
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
                    TS, TSextra, nbins=c(50, 10), rate=0, 
                    Range=c(-Inf, Inf), B=c(5000,1000),  minexpcount=5.0,  
                    ChiUsePhat = FALSE, maxProcessor, 
                    SuppressMessages =FALSE, seed=NULL, doMethods, qnull) {
 
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
   # qnull can be given as an argument or (as in earlier versions) in TSextra
   if(!missing(qnull)) {
     if(missing(TSextra)) TSextra <- list()
     TSextra$qnull <- qnull
   }
   TSextra=makeTSextra(TSextra, pnull, phat, w, Continuous)
   WithWeights = TRUE
   if(length(formals(w))==1) {
     if(w(x[1])==-99) WithWeights = FALSE
   }
   # adjust number of bins to account for parameter estimation
   if(abs(phat(x)[1]+99)>0.001) nbins=nbins+length(phat(x))
   tmp <- maketypeTS(TS, Continuous, WithWeights)
   typeTS <- tmp$typeTS
   TS <- tmp$TS
   useSingleProcessor <- tmp$useSingleProcessor
   if(!is.null(seed)) {
     if(length(seed)!=1L || is.na(seed) || !is.finite(seed))
       stop("seed must be NULL or a single finite number", call.=FALSE)
     
     set.seed(seed)
   }
   TS_data=calcTS(dta, TS, typeTS, TSextra)  
   if(is.null(names(TS_data))) {
     message("result of TS has to be a named vector")
     return(NULL)
   }
   # Choose number of processors after dta, TS and typeTS are known.
   maxProcessor <- makemaxProcessor(maxProcessor,
             dta=dta, TS=TS, typeTS=typeTS, TSextra=TSextra, B=B,
             SuppressMessages=SuppressMessages,
             useSingleProcessor=useSingleProcessor
   )
   B[1]=ceiling(B[1]/maxProcessor)*maxProcessor
   if(length(x)>10000 && maxProcessor==1)
     message("Consider using parallel processing with maxProcessor= (your number of cores)")
   
   if(maxProcessor==1) {
      A=gof_adj_C1(dta, rnull, vals, TS, typeTS, TSextra, B[1])
      pvals=gof_adj_C2(dta, rnull, vals, TS, typeTS, TSextra, A, B[2])
   }    
   else {
     cl = parallel::makeCluster(maxProcessor)
     on.exit(parallel::stopCluster(cl), add=TRUE)
     if(!is.null(seed))
       parallel::clusterSetRNGStream(cl, iseed=seed)
     z=parallel::clusterCall(cl, gof_adj_C1, 
                dta, rnull, vals, TS, typeTS, TSextra, 
                             B=round(B[1]/maxProcessor)
     )
     A=z[[1]]
     for(i in 2:maxProcessor) A=rbind(A, z[[i]])
     z=parallel::clusterCall(cl, gof_adj_C2, 
                          dta, rnull, vals, TS, typeTS, TSextra, 
                          A, B=round(B[2]/maxProcessor)
     )
# row 1 of each result is the observed data, keep it only once
     pvals=z[[1]]
     for(i in 2:maxProcessor) pvals=rbind(pvals,z[[i]][-1, , drop=FALSE])
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
      pvals=cbind(pvals, chipvals)
   }
   # without qnull the Wasserstein test can not be done for continuous data
   if(Continuous && TSextra$Noqnull) {
     pvals = pvals[, colnames(pvals)!="Wassp1", drop=FALSE]
     if(!missing(doMethods)) {
       if("Wassp1" %in% doMethods && !SuppressMessages)
         message("Wassp1 needs qnull and is not run")
       doMethods = doMethods[doMethods!="Wassp1"]
     }
   }
   if(missing(doMethods)) doMethods=colnames(pvals)
   if(test_methods(doMethods, Continuous, WithWeights)) return(NULL)
   pvals = round(pvals[, doMethods, drop=FALSE],4)
   minp_x=min(pvals[1, ])
   minp_sim=apply(pvals[-1, ,drop=FALSE], 1, min)
# adjusted p value = P(smallest p value <= observed smallest p value | H0)
   minp_adj <- round(mean(minp_sim<=minp_x), 4)
   message("p values of individual tests:")
   for(i in 1:ncol(pvals)) message(paste(doMethods[i],": ", pvals[1,i]))
   message(paste0("adjusted p value of combined tests: ", minp_adj))
   invisible(list(p.values=stats::setNames(pvals[1, ], colnames(pvals)), adjusted.p.value=minp_adj))

}

#' This function performs a number of gof tests for discrete data and finds the adjusted p value
#' @param  x data set (the counts)
#' @param  pnull  cumulative distribution function under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  vals a vector of values of discrete random variables 
#' @param  phat =function(x) -99, function to estimate parameters from the data, or -99 if no parameters aare estimated
#' @param  TS =NA, user supplied function to find test statistics
#' @param  TSextra =NA, list passed to TS, if desired
#' @param  nbins =c(50, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  B   =c(5000, 1000)  number of simulation runs for p values and for adjusted p value
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  doMethods Methods to include in tests
#' @return A numeric matrix of test statistics and p.values

gof_test_disc_adj <- function(x, pnull, rnull, vals, phat=function(x) -99, 
                          TS, TSextra=NA,  nbins=c(50, 10), rate=0, 
                          B=c(5000, 1000), minexpcount=5.0, ChiUsePhat=TRUE,
                          doMethods=c("Wassp1", "W", "AD", "s-P")) {
  if(length(B)==1) B=c(B, B)
  if(any(is.na(TSextra))) TSextra = list(p=phat(x))
  else TSextra = c(TSextra, p=phat)
  if(missing(TS)) { # use built-in tests
    typeTS = 5
    TS = TS_disc
    TS_data = TS(x, pnull, phat(x), vals)
  }  
  else {
    # can't do parallel processing if TS written in C/C++
    if(substr(deparse(TS)[2], 1, 5)==".Call") {
      message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
      maxProcessor=1
    }
    if(length(formals(TS))==4) {
      typeTS=5
      TS_data = TS(x, pnull, phat(x), vals)
    }
    if(length(formals(TS))==5) {
      typeTS=6
      TS_data = TS(x, pnull, phat(x), vals, TSextra)
    }
    if(length(formals(TS))>6) {
      message("For discrete data TS should have either 4 or 5 arguments")
      return(NULL)
    }
    if(is.null(names(TS_data))) {
      message("result of TS has to be a named vector")
      return(NULL)
    } 
  }
  NoEstimation=FALSE
  if(length(formals(pnull))==0) NoEstimation=TRUE
  if(typeTS>5) doMethods=names(TS_data)
  
  num_tests=length(TS_data)
  A=matrix(0, B[1], num_tests)
  p=phat(x)
  for(i in 1:B[1]) {
     if(NoEstimation) xsim=rnull()
     else xsim=rnull(p)
     TS_sim=calcTS(list(x=xsim, vals=vals), pnull, phat(xsim), TS, typeTS, TSextra)      
     A[i, ]=TS_sim
  }
  if(typeTS==5) {
      pvals=matrix(0, B[2]+1, num_tests+4)
      colnames(pvals)=c(names(TS_data), "l-P", "s-P", "l-L", "s-L")      
  }    
  else {
    pvals=matrix(0, B[2]+1, num_tests)
    colnames(pvals)=names(TS_data)
  }  
  for(i in 1:(B[2]+1)) {
     if(i==1) xsim=x
     else {
          if(NoEstimation) xsim=rnull()
          else xsim=rnull(p) 
     }      
     TS_sim=calcTS(list(x=xsim, vals=vals), pnull, phat(xsim), TS, typeTS, TSextra)
     for(j in 1:num_tests) 
        pvals[i, j]=pvals[i, j]+sum(TS_sim[j]<A[,j])/B[1]
     if(typeTS==5) {
        pvals[i, num_tests+1:4]=round(chi_test_disc(xsim, pnull, 
                         phat=phat, nbins=nbins, rate=0, 
                         minexpcount=minexpcount, 
                         ChiUsePhat=ChiUsePhat)[,2],4)
       }                  
  }
  pvals=pvals[, doMethods, drop=FALSE]
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

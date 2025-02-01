#' This function performs a number of gof tests for continuous data and finds the adjusted p value
#' @param  x data set
#' @param  pnull  cdf under the null hypothesis
#' @param  rnull  routine to generate data under the null hypothesis
#' @param  w (Optional) function to calculate weights, returns -99 if no weights
#' @param  phat =function(x) -99, function to estimate parameters from the data, or -99 if no parameters aare estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  TSextra =NA, list passed to TS, if desired
#' @param  nbins =c(50, 10) number of bins for chi-square tests
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any, for chi-square tests
#' @param  B   =c(5000,1000)  number of simulation runs for p values and for p value distribution
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat =TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @param  doMethods Methods to include in tests
#' @return None


gof_test_cont_adj=function(x, pnull, rnull, w=function(x) -99, phat=function(x) 0, 
       TS, TSextra=NA, nbins=c(50,10), rate=0, Range=c(-Inf,Inf), 
       B=c(5000, 1000), minexpcount=5, ChiUsePhat=TRUE, 
       doMethods=c("W", "ZC", "AD", "ES-s-P")) {
  if(length(B)==1) B=c(B, B)
 # Are weights present?
  WithWeights = TRUE
  if(length(formals(w))==1 & w(x[1])==-99) WithWeights = FALSE
  if(any(is.na(TSextra))) TSextra = list(p=phat(x))
  else TSextra = c(TSextra, p=phat)
  Noqnull = FALSE
  if( !("qnull" %in% names(TSextra)) ) {
    Noqnull = TRUE
    qnull=function(x, p=0) rep(-99,length(x))
    TSextra = c(TSextra, qnull=qnull)
  }  
  else qnull = TSextra$qnull      
  if(missing(TS)) {
     if(!WithWeights) { #data is not weighted
       typeTS=1
       TS = TS_cont
       TS_data = TS(x, pnull, phat(x), function(x) abs(x)/max(x))
     }
     else {
       typeTS=2
       TS = TSw_cont
       TS_data = TS(x, pnull, phat(x), w(x))
       doMethods = names(TS_data)
     }
  }   
  else {
    if(length(formals(TS))==3) {
       typeTS=3
       TS_data = TS(x, pnull, phat(x))
    }
    if(length(formals(TS))==4) {
      typeTS=4
      TS_data = TS(x, pnull, phat(x), TSextra)
    }
    if(length(formals(TS))>4) {
      message("TS should have either 3 or 4 arguments")
      return(NULL)
    }
    if(is.null(names(TS_data))) {
      message("result of TS has to be a named vector")
      return(NULL)
    }
  }
  p=phat(x)
  psim=p
  NoEstimation=FALSE
  if(length(formals(pnull))==1) NoEstimation=TRUE
  if(WithWeights) {
    if(NoEstimation) wx=w(x)   
    else wx=w(x,p)
  }
  if(typeTS==1) TS_data=TS(x, pnull, p, qnull);
  if(typeTS==2) TS_data=TS(x, pnull, p, wx);  
  if(typeTS==3) TS_data=TS(x, pnull, p);
  if(typeTS==4) TS_data=TS(x, pnull, p, TSextra);
  if(typeTS>2) doMethods=names(TS_data)
  num_tests=length(TS_data)
  A=matrix(0, B[1], num_tests)
  for(i in 1:B[1]) {
     if(NoEstimation) xsim=rnull()
     else {xsim=rnull(p);psim=phat(xsim)}
     if(WithWeights) {
       if(NoEstimation) wx=w(xsim)   
       else wx=w(xsim, psim)
     }
     if(typeTS==1) TS_sim=TS(xsim, pnull, psim, qnull);
     if(typeTS==2) TS_sim=TS(xsim, pnull, psim, wx);  
     if(typeTS==3) TS_sim=TS(xsim, pnull, psim);
     if(typeTS==4) TS_sim=TS(xsim, pnull, psim, TSextra);    
     A[i, ]=TS_sim
  }
  if(typeTS<=2) {
      pvals=matrix(0, B[2]+1, num_tests+8)
      colnames(pvals)=c(names(TS_data),
                    "ES-l-P", "ES-s-P", "EP-l-P", "EP-s-P",
                    "ES-l-L", "ES-s-L", "EP-l-L", "EP-s-L")      
  }    
  else {
    pvals=matrix(0, B[2]+1, num_tests)
    colnames(pvals)=names(TS_data)
  }  
  for(i in 1:(B[2]+1)) {
     if(i==1) {xsim=x;psim=p}
     else {
          if(NoEstimation) xsim=rnull()
          else {xsim=rnull(p);psim=phat(xsim)}
     }      
     if(NoEstimation) {
         if(typeTS==2) wx=w(xsim)   
     }    
     else {
         if(typeTS==2) wx=w(xsim, psim)
     }
     if(typeTS==1) TS_sim=TS(xsim, pnull, psim, qnull);
     if(typeTS==2) TS_sim=TS(xsim, pnull, psim, wx);  
     if(typeTS==3) TS_sim=TS(xsim, pnull, psim);
     if(typeTS==4) TS_sim=TS(xsim, pnull, psim, TSextra);
     for(j in 1:num_tests) 
        pvals[i, j]=pvals[i, j]+sum(TS_sim[j]<A[,j])/B[1]
     if(typeTS<=2) {
        Range[1]=ifelse(is.infinite(Range[1]),-99999, Range[1])
        Range[2]=ifelse(is.infinite(Range[2]),99999, Range[2])
        pvals[i, num_tests+1:8]=round(chi_test_cont(xsim, pnull, w=w,
                phat=phat, qnull=ifelse(Noqnull, NA, qnull), rate=rate, nbins=nbins,
                         Range=Range, minexpcount=minexpcount, 
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

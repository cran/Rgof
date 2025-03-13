#'  This function runs the case studies included in the package
#' @param TS routine to calculate test statistic(s) or p value(s).
#' @param study either the name of the study, or its number. If missing all the studies are run.
#' @param TSextra =list(aaa=1), list passed to TS.
#' @param With.p.value =FALSE does user supplied routine return p values?
#' @param BasicComparison =TRUE if true compares tests on one default value of parameter of the alternative distribution.
#' @param nsample = 500, desired sample size.
#' @param alpha =0.05  type I error
#' @param param_alt (list of) values of parameter under the alternative hypothesis. If missing included values are used.
#' @param maxProcessor number of cores to use for parallel programming
#' @param B = 1000 number of simulation runs
#' @return A (list of ) matrices of p.values
#' @examples
#' # New test is a simple chi-square test: 
#' chitest=function(x, pnull, param, TSextra) {
#'     nbins=TSextra$nbins
#'     bins=quantile(x, (0:nbins)/nbins)
#'     O=hist(x, bins, plot=FALSE)$counts
#'     if(param[1]!=-99) { #with parameter estimation
#'         E=length(x)*diff(pnull(bins, param))
#'         chi=sum((O-E)^2/E)
#'         pval=1-pchisq(chi, nbins-1-length(param))
#'     }
#'     else {
#'       E=length(x)*diff(pnull(bins))
#'       chi=sum((O-E)^2/E)
#'       pval=1-pchisq(chi,nbins-1)
#'     }  
#'     out=ifelse(TSextra$statistic, chi, pval)
#'     names(out)="ChiSquare"
#'     out
#' }
#' TSextra=list(nbins=10, statistic=FALSE) # Use 10 bins, test routine returns p-value
#' run.studies(chitest, TSextra=TSextra, With.p.value=TRUE, maxProcessor=1, B=200)
#' @export

run.studies <- function(TS, study, TSextra=list(aaa=1), With.p.value=FALSE, BasicComparison=TRUE, 
                nsample=500, alpha=0.05, param_alt, maxProcessor, B=1000) {
  
  B=B[1] 
  if(!is.function(TS)) {
      if(missing(TS) || !is.logical(TS)) {
         message("TS should either be a function, or TRUE/FALSE (for continuous or discrete) to run the included tests") 
         return(NULL)
      }   
      Continuous=TS
  }    
  else {
    Continuous=ifelse("vals" %in% names(formals(TS)), FALSE, TRUE)  
    WithTSextra=ifelse("aaa"%in%names(TSextra), FALSE, TRUE) 
    if(Continuous) {
      if(any(names(unlist(formals(TS))[1:3])!=c("x", "pnull", "param"))) {
          message("for continuous data the function TS should have arguments x, pnull, param and TSextra (optional)")
          return(NULL)
      }    
    }
    else {
      if(any(names(unlist(formals(TS))[1:4])!=c("x", "pnull", "param", "vals"))) {
         message("for discrete data the function TS should have arguments x, pnull, param, vals and TSextra (optional)")
         return(NULL)
      }    
    }
  }   
  if(missing(maxProcessor)) maxProcessor=parallel::detectCores()-1
  I80cont=c(12,11,8,11,13,7,9,9,24,8,12,13,7,13,11,19,8,14,17,13)
  I80disc=c(12,11,9,11,14,7,9,9,21,9,11,13,8,15,8,23,8,14,15,12)
  if(Continuous) I80=I80cont
  else I80=I80disc
  NewParams=ifelse(alpha==0.05, FALSE, TRUE)
  if(!missing(param_alt)) {
      NewParams=TRUE
      message("For new parameter values under the alternative or alpha!=0.05 power values 
        will also be calculated for included tests")
      if(!missing(study) && length(study)>1) {
        if(!is.list(param_alt)) {
           message("param_alt has to be a list with the same length as study")
           return(NULL)
        }
      }
      else if(!is.list(param_alt)) param_alt=list(param_alt)
  } 
  if(NewParams) BasicComparison=FALSE 
  list.of.studies=c(
  "uniform.linear",
  "uniform.quadratic",
  "uniform.bump",
  "uniform.sine",
  "beta22.betaaa",
  "beta22.beta2a",
  "normal.shift",
  "normal.stretch",
  "normal.t",
  "normal.outlier1",
  "normal.outlier2",
  "exponential.gamma",
  "exponential.weibull",
  "exponential.bump",
  "trunc.exponential.linear",
  "normal.t.est",
  "exponential.weibull.est",
  "trunc.exponential.linear.est",
  "exponential.gamma.est",
  "normal.cauchy.est"
  )
  if(Continuous) list.of.studies=paste0(list.of.studies,".cont")
  else list.of.studies=paste0(list.of.studies,".disc")
  if(missing(study)) study=1:20
  else BasicComparison=FALSE
  if(is.numeric(study)) study=list.of.studies[study]
  else {
    if(Continuous) study=paste0(study,".cont")
    else study=paste0(study,".disc")
  }
  out=as.list(seq_along(study))
  names(out)=study
  for(i in seq_along(study)) {
    message(paste("Running case", study[i], "..."))
    WithEstimation=FALSE
    if(endsWith(substring(study[i],1,nchar(study[i])-5),"est")) WithEstimation=TRUE
    pwrold=Rgof::power_studies_results[[study[i]]]
    tmp=case.studies(study[i], nsample)
    TSextra$pnull=tmp$pnull
    if("phat"%in%names(tmp)) phat=tmp$phat
    else phat=function(x) -99
    if(BasicComparison) tmp$param_alt=tmp$param_alt[I80[i]]
    if(NewParams || !is.function(TS)) {
       if(!missing(param_alt)) tmp$param_alt=param_alt[[i]]
       if(Continuous) {
            if(WithEstimation) pwrold=Rgof::gof_power(tmp$pnull, NA, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, phat=tmp$phat, maxProcessor=maxProcessor,B=B)
            else pwrold=Rgof::gof_power(tmp$pnull, NA, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, maxProcessor=maxProcessor,B=B)            
       }
       else {
            if(WithEstimation) pwrold=Rgof::gof_power(tmp$pnull, tmp$vals, tmp$rnull, 
                    tmp$ralt, tmp$param_alt, phat=tmp$phat, maxProcessor=maxProcessor,B=B)
            else pwrold=Rgof::gof_power(tmp$pnull, tmp$vals, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, maxProcessor=maxProcessor,B=B) 
       }
       if(!is.matrix(pwrold)) {
          pwrold=rbind(pwrold)
          rownames(pwrold)=tmp$param_alt
       }   
    }
    if(!is.function(TS)) {out[[i]]=pwrold;next} 
    if(With.p.value) {
        if(Continuous) {
           pwr=power_newtest(TS, NA, tmp$pnull,
               tmp$ralt, tmp$param_alt, phat, TSextra, alpha, B[1])     
        }
        else {
            pwr=power_newtest(TS, tmp$vals, tmp$pnull,
               tmp$ralt, tmp$param_alt, phat, TSextra, alpha, B[1])     
        } 
    }
    else {
        if(Continuous && WithEstimation && WithTSextra) 
               pwr=Rgof::gof_power(tmp$pnull, NA, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, phat=tmp$phat, TS=TS, 
                          TSextra=TSextra, maxProcessor=maxProcessor, B=B)
        if(Continuous && !WithEstimation && WithTSextra)
               pwr=Rgof::gof_power(tmp$pnull, NA, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, TS=TS, TSextra=TSextra, 
                                maxProcessor=maxProcessor, B=B)
        if(!Continuous && WithEstimation && WithTSextra)
              pwr=Rgof::gof_power(tmp$pnull, tmp$vals, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, phat=tmp$phat, TS=TS, 
                    TSextra=TSextra, maxProcessor=maxProcessor, B=B)
        if(!Continuous && !WithEstimation && WithTSextra)
              pwr=Rgof::gof_power(tmp$pnull, tmp$vals, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, TS=TS, TSextra=TSextra, 
                                maxProcessor=maxProcessor, B=B)
        if(Continuous && WithEstimation && !WithTSextra) 
               pwr=Rgof::gof_power(tmp$pnull, NA, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, phat=tmp$phat, TS=TS, 
                                maxProcessor=maxProcessor, B=B)
        if(Continuous && !WithEstimation && !WithTSextra)
               pwr=Rgof::gof_power(tmp$pnull, NA, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, TS=TS, 
                                maxProcessor=maxProcessor, B=B)      
        if(!Continuous && WithEstimation && !WithTSextra)
              pwr=Rgof::gof_power(tmp$pnull, tmp$vals, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, phat=tmp$phat, TS=TS, 
                                maxProcessor=maxProcessor, B=B)
        if(!Continuous && !WithEstimation && !WithTSextra)
              pwr=Rgof::gof_power(tmp$pnull, tmp$vals, tmp$rnull, 
                   tmp$ralt, tmp$param_alt, TS=TS, 
                                maxProcessor=maxProcessor, B=B)   
        if(length(tmp$param_alt)==1) {       
            dta=tmp$ralt(tmp$param_alt[1])
            if(Continuous && !WithTSextra) nm=TS(dta, tmp$pnull, phat(dta)) 
            if(Continuous && WithTSextra) nm=TS(dta, tmp$pnull, phat(dta), TSextra) 
            if(!Continuous && !WithTSextra) nm=TS(dta, tmp$pnull, phat(dta), tmp$vals) 
            if(!Continuous && WithTSextra) nm=TS(dta, tmp$pnull, phat(dta), tmp$vals, TSextra)
            pwr=matrix(pwr, 1, length(pwr))     
            colnames(pwr)=names(nm)
            rownames(pwr)=tmp$param_alt
        }        
    }   
    if(BasicComparison) out[[i]]=cbind(pwr, pwrold[I80[i], , drop=FALSE])
    else out[[i]]=cbind(pwr[, , drop=FALSE], pwrold[, , drop=FALSE])
  } 
  if(BasicComparison) {
     A=matrix(0, 20, ncol(out[[1]]))
     for(i in 1:20) A[i, ]= out[[i]][1, ]
     rownames(A) = list.of.studies
     colnames(A) = colnames(out[[1]])
     a1=apply(A, 1, rank)
     message("Average rank of the tests:")
     print(sort(apply(a1,1,mean)))
     return(A)
  }
  if(length(out)==1) return(out[[1]])
  out
}

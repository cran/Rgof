#' Power estimation of goodness-of-fit tests.
#' 
#' Find the power of various goodness-of-fit tests.
#' 
#' For details on the usage of this routine consult the vignette with vignette("Rgof","Rgof")
#' 
#' @param  pnull function to find cdf under  null hypothesis
#' @param  vals =NA values of discrete random variable, or NA
#' @param  rnull function to generate data under  null hypothesis
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  w (Optional) function to calculate weights, returns -99 if no weights
#' @param  phat =function(x) -99 function to estimate parameters from the data, or -99
#' @param  TS user supplied function to find test statistics
#' @param  TSextra list provided to TS (optional)
#' @param  With.p.value =FALSE does user supplied routine return p values?
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  Range  =c(-Inf, Inf) limits of possible observations, if any
#' @param  B =1000 number of simulation runs
#' @param  nbins =c(50,10), number of bins for chi square tests.
#' @param  rate =0 rate of Poisson if sample size is random, 0 if sample size is fixed
#' @param  maxProcessor maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @param  minexpcount =5 minimal expected bin count required
#' @param  ChiUsePhat = TRUE, if TRUE param is estimated parameter, otherwise minimum chi square method is used.
#' @return A numeric matrix of power values.
#' @export 
#' @examples
#' # Power of tests when null hypothesis specifies the standard normal distribution but 
#' # true data comes from a normal distribution with mean different from 0.
#' pnull = function(x) pnorm(x)
#' rnull = function()  rnorm(50)
#' ralt = function(mu)  rnorm(50, mu)
#' TSextra = list(qnull=function(x) qnorm(x))
#' gof_power(pnull, NA, rnull, ralt, c(0.25, 0.5), TSextra=TSextra, B=200)
#' # Power of tests when null hypothesis specifies normal distribution and 
#' # mean and standard deviation are estimated from the data. 
#' # true data comes from a normal distribution with mean different from 0.
#' pnull = function(x, p=c(0, 1)) pnorm(x, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' rnull = function(p=c(0, 1))  rnorm(50, p[1], ifelse(p[2]>0.001, p[2], 0.001))
#' ralt = function(mu)  rnorm(50, mu)
#' phat = function(x) c(mean(x), sd(x))
#' TSextra = list(qnull = function(x, p=c(0, 1)) qnorm(x, p[1],  
#'                ifelse(p[2]>0.001, p[2], 0.001))) 
#' pwr=gof_power(pnull, NA, rnull, ralt, c(0, 1), phat=phat, TSextra=TSextra, B=200)
#' pwr
#' #' Compare power of a new test based on variants of the Cramer-vonMises
#' #' criterion to the methods included in the package: 
#' newTS = function(x, pnull, param) {
#'    Fx=sort(pnull(x, param))
#'    n=length(x)
#'    out = c(sum(abs( (2*1:n-1)/2/n-Fx )), sum(sqrt(abs( (2*1:n-1)/2/n-Fx ))))
#'    names(out) = c("CvM alt 1","CvM alt 2")
#'    out
#' }
#' #' Compare power to Lilliefors KS test, which finds its own p value:
#' LLtest=function(x, pnull, param) {
#'   out=nortest::lillie.test(x)$p.value
#'   names(out)="KS - Lilliefors"
#'   out
#' }
#' cbind(gof_power(pnull, NA, rnull, ralt, c(0, 1), TS=LLtest, phat=phat, 
#'        With.p.value=TRUE, TSextra=TSextra, B=200), pwr)
#' # Power of tests when null hypothesis specifies Poisson rv with rate 100 and 
#' # true rate is 100.5
#' vals = 0:250
#' pnull = function() ppois(0:250, 100)
#' rnull =function () table(c(0:250, rpois(1000, 100)))-1
#' ralt =function (p) table(c(0:250, rpois(1000, p)))-1
#' gof_power(pnull, vals, rnull, ralt, param_alt=100.5,  B=200)
#' # Power of tests when null hypothesis specifies a Binomial n=10 distribution 
#' # with the success probability estimated
#' vals = 0:10
#' pnull=function(p) pbinom(0:10, 10, ifelse(0<p&p<1, p, 0.001))
#' rnull=function(p) table(c(0:10, rbinom(1000, 10, ifelse(0<p&p<1, p, 0.001))))-1
#' ralt=function(p) table(c(0:10, rbinom(1000, 10, p)))-1
#' phat=function(x) mean(rep(0:10,x))/10
#' gof_power(pnull, vals, rnull, ralt, c(0.5, 0.6), phat=phat, B=200)
#'
gof_power=function(pnull, vals=NA, rnull, ralt, param_alt, 
        w=function(x) -99, phat=function(x) -99, TS, TSextra, 
        With.p.value=FALSE, 
        alpha=0.05, Range  =c(-Inf, Inf), B=1000,nbins=c(50,10), 
        rate=0, maxProcessor, minexpcount=5.0, ChiUsePhat=TRUE) {

  fff=nortest::lillie.test # avoid issues with CRAN, just ignore!
  NewTest=TRUE
  if(missing(TS)) NewTest=FALSE
  dta = ralt(param_alt[1]) # get an example data set
  x = dta
  Continuous=ifelse(any(is.na(vals)), TRUE, FALSE)
  if(Continuous) {
    dta=list(x=x)
    check.functions(pnull, rnull, phat, x=x)
  }  
  else {
    dta=list(x=x, vals=vals)
    check.functions(pnull, rnull, phat, vals, x)
  }

  if(missing(TSextra)) TSextra=list(pnull=pnull, phat=phat, 
                            w=w, Continuous=Continuous)
  else TSextra = c(TSextra, pnull=pnull, phat=phat, 
                   w=w, Continuous=Continuous)
  Noqnull=FALSE
  if(!("qnull" %in% names(TSextra))) {
    Noqnull=TRUE
    TSextra=c(TSextra, qnull=function(x) -99)
  } 
  WithWeights = TRUE
  if(length(formals(w))==1) {
    if(w(x[1])==-99) WithWeights = FALSE
  }
  # adjust number of bins to account for parameter estimation
  if(abs(phat(x)[1]+99)>0.001) nbins=nbins+length(phat(x))
  if(any(is.na(vals))) check.functions(pnull, rnull, phat, x=x)
  else  check.functions(pnull, rnull, phat, vals, x)
  if(missing(TS)) {
    if(Continuous) {
      if(!WithWeights) { #data is not weighted
        typeTS=1
        TS = TS_cont
      }
      else {
        typeTS=2
        TS = TSw_cont
      }
    }
    else {
      typeTS = 5
      TS = TS_disc
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
  if(missing(maxProcessor)) 
    maxProcessor=parallel::detectCores(logical = FALSE)-1
    if(With.p.value) maxProcessor=1
  if(maxProcessor>1) {
    tm=timecheck(dta, TS, typeTS, TSextra)
    if(tm*length(param_alt)*B<20) {
      maxProcessor=1
      message("maxProcessor set to 1 for faster computation")
    }
    else message(paste("Using ",maxProcessor," cores.."))
  }
  if(With.p.value) {
    if(Continuous) {
      pwr=power_newtest(TS, NA, pnull, ralt, param_alt, TSextra$phat, TSextra, alpha, B)     
    }
    else {
      pwr=power_newtest(TS, vals, pnull, ralt, param_alt, TSextra$phat, TSextra, alpha, B)     
    } 
  }
  else {
    if(maxProcessor==1) {
        tmp=gof_power_C(rnull, vals, ralt, param_alt, TS, typeTS, TSextra, B)
        Data=tmp$Data
        Sim=tmp$Sim
    }
    else {
        cl <- parallel::makeCluster(maxProcessor)
        z=parallel::clusterCall(cl, gof_power_C, 
                rnull, vals, ralt, param_alt,  TS, typeTS, TSextra, 
                B=round(B/maxProcessor))
        parallel::stopCluster(cl)
        Sim=z[[1]][["Sim"]]
        Data=z[[1]][["Data"]]
        for(i in 2:maxProcessor) {
          Sim=rbind(Sim,z[[i]][["Sim"]])
          Data=rbind(Data,z[[i]][["Data"]])
        }  
    }  
    pwr=matrix(0, length(param_alt), length(TS_data))
    colnames(pwr)=names(TS_data)
    rownames(pwr)=param_alt
    crtval=apply(Data, 2, quantile, prob=1-alpha, na.rm=TRUE)
    for(i in seq_along(param_alt)) {
      tmpS=Sim[Sim[,1]==param_alt[i], -1, drop=FALSE]
      for(j in seq_along(crtval)) 
        pwr[i, j]=sum(tmpS[ ,j]>crtval[j])/nrow(tmpS)
    }
  }
  # Do chi square tests if built-in TS is used. Don't run chi square if weights are present.  
  chipwr=NULL
  if(typeTS==1) { #Run chi square tests
    if(is.infinite(Range[1])) Range[1]=-99999
    if(is.infinite(Range[2])) Range[2]=99999  
    chipwr = chi_power_cont(pnull=pnull, 
                            ralt = ralt, 
                            param_alt = param_alt,                             
                            qnull = ifelse(Noqnull, NA, TSextra$qnull), 
                            phat = phat, 
                            w = w,
                            alpha = alpha, 
                            Range = Range, 
                            B= B, 
                            nbins = nbins, 
                            rate = rate, 
                            minexpcount = minexpcount,
                            ChiUsePhat=ChiUsePhat) 
  }
  if(typeTS==5 & (!NewTest)) { #Run chi square tests
    chipwr = chi_power_disc(pnull, ralt, param_alt, 
                            phat, alpha , B, 
                            nbins, rate, minexpcount,
                            ChiUsePhat)[,1:2, drop=FALSE]
  }
  if(typeTS==1 | typeTS==5) pwr = cbind(pwr, chipwr)
  if(is.matrix(pwr) & nrow(pwr)==1) pwr=pwr[1, ]
  round(pwr, 3)
}

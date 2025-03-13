#' Find the power of various gof tests for discrete data.
#' @param  pnull function to find cdf under  null hypothesis
#' @param  rnull function to generate data under  null hypothesis
#' @param  vals  values of discrete distribution
#' @param  ralt function to generate data under  alternative hypothesis
#' @param  param_alt  vector of parameter values for distribution under alternative hypothesis
#' @param  phat =function(x) -99, function to estimate parameters from the data, or -99 if no parameters are estimated
#' @param  TS user supplied function to find test statistics, if any
#' @param  typeTS format of TS routine
#' @param  TSextra list provided to TS
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B =1000 number of simulation runs
#' @param  maxProcessor maximum of number of processors to use, 1 if no parallel processing is needed or number of cores-1 if missing
#' @return A numeric matrix of power values

power_disc_R=function(pnull, rnull, vals, ralt, param_alt, 
        phat=function(x) -99, TS, typeTS, TSextra, 
        alpha=0.05,  B=1000,  maxProcessor) {
  
  x = ralt(param_alt[1])  
  WithEstimation=0
  if(abs(phat(x)[1]+99)>0.001) WithEstimation=1
  TS_data=calcTS(list(x=x, vals=vals), pnull, phat(x), TS, typeTS, TSextra)
  nummethods = length(TS_data)
  methods = names(TS_data)

  if(maxProcessor==1) { # no parallel processing 
    tmp=power_disc(pnull, rnull, vals, ralt, param_alt, phat, 
                   TS, typeTS, TSextra, B)
    Data=tmp$Data
    Sim=tmp$Sim
  }  
  else { 
    cl <- parallel::makeCluster(maxProcessor)
    z=parallel::clusterCall(cl, power_disc, 
                    pnull, rnull, vals, ralt, param_alt, phat, 
                    TS, typeTS, TSextra, B=round(B/maxProcessor))
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
  round(pwr, 3)
}


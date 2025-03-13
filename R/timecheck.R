#' estimate run time function
#' @param  x data set
#' @param  pnull function to find cdf under  null hypothesis
#' @param  phatx parameter estimates
#' @param  wx vector of wights
#' @param  TS test statistic
#' @param  typeTS format of TS
#' @param  TSextra additional info TS
#' @return Mean computation time
#' @export
timecheck=function(x, pnull, phatx, wx, TS, typeTS, TSextra) {
  if(typeTS==1) f=function() TS(x, pnull, phatx, TSextra$qnull);
  if(typeTS==2) f=function() TS(x, pnull, phatx, wx);
  if(typeTS==3) f=function() TS(x, pnull, phatx);
  if(typeTS==4) f=function() TS(x, pnull, phatx, TSextra);
  a=microbenchmark::microbenchmark(f(), 
              unit="seconds", times=10)
  as.numeric(summary(a)["mean"])
}

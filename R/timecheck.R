#' estimate run time function
#' @param  dta data set
#' @param  TS test statistic
#' @param  typeTS format of TS
#' @param  TSextra additional info TS
#' @return Mean computation time
#' @export
timecheck=function(dta, TS, typeTS, TSextra) {
  f=function() calcTS(dta, TS, typeTS, TSextra)
  a=microbenchmark::microbenchmark(f(), 
              unit="seconds", times=10)
  as.numeric(summary(a)["mean"])
}

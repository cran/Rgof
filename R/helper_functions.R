#' makeTSextra
#' 
#' create the TSextra list required by many routines.
#' 
#' @param TSextra a user-supplied list
#' @param pnull cdf
#' @param phat estimation routine
#' @param w weight function
#' @param Continuous is data continuous?
#' @return a list
#' @keywords internal
#' @export
makeTSextra=function(TSextra, pnull, phat, w, Continuous) {
  if(missing(TSextra)) TSextra <- list()
  TSextra$pnull <- pnull
  TSextra$phat <- phat
  TSextra$w <- w
  TSextra$Continuous <- Continuous
  Noqnull=FALSE
  if(!("qnull" %in% names(TSextra))) {
    Noqnull=TRUE
    TSextra=c(TSextra, qnull=function(x) -99)
  }
  TSextra$Noqnull=Noqnull
  TSextra
}
#' maketypeTS
#' 
#' find typeTS and TS
#' 
#' @param TS a function
#' @param Continuous is data continuous?
#' @param WithWeights with weights?
#' @return a list
#' @keywords internal
#' @export
maketypeTS=function(TS, Continuous, WithWeights) { 
  useSingleProcessor <- FALSE
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
      useSingleProcessor <- TRUE
    }
    nargs <- length(formals(TS))
    if(Continuous) {
      if(!(nargs %in% c(3L,4L))) {
        stop("TS for continuous data should have either 3 or 4 arguments",
           call.=FALSE)
      }
      typeTS=nargs
    }  
    else {
      if(!(nargs %in% c(4L,5L))) {
        stop("TS for discrete data should have either 4 or 5 arguments",
           call.=FALSE)
      }
      typeTS=nargs+1
    }
  }  
  list(TS=TS, typeTS=typeTS, 
       useSingleProcessor=useSingleProcessor)
}

#' validate_integer_scalar
#'
#' Validate a positive integer-valued scalar.
#'
#' @param x value to validate
#' @param name argument name used in error messages
#' @param lower smallest permitted value
#' @return invisibly returns the validated integer
#' @keywords internal
#' @noRd
validate_integer_scalar <- function(x, name, lower=1L) {
  if(!is.numeric(x) || length(x)!=1L || is.na(x) || !is.finite(x) ||
     x != floor(x) || x < lower)
    stop(name, " must be a single integer >= ", lower, call.=FALSE)
  invisible(as.integer(x))
}

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
#' @noRd
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
#' @noRd
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

#' makemaxProcessor
#'
#' Choose the number of processors for a simulation.
#'
#' @param maxProcessor maximum number of processors requested
#' @param dta data set used for a timing check
#' @param TS test-statistic routine
#' @param typeTS test-statistic interface code
#' @param TSextra additional information passed to TS
#' @param B number of simulation runs
#' @param SuppressMessages suppress informative messages?
#' @param useSingleProcessor force a single processor, e.g. for a C++ custom TS
#' @return an integer giving the number of processors to use
#' @keywords internal
#' @noRd
makemaxProcessor <- function(maxProcessor, dta, TS, typeTS, TSextra, B,
                             SuppressMessages=FALSE,
                             useSingleProcessor=FALSE) {
  if(useSingleProcessor) return(1L)
  
  if(!missing(maxProcessor)) {
    validate_integer_scalar(maxProcessor, "maxProcessor", 1L)
    return(as.integer(maxProcessor))
  }
  
  m <- parallel::detectCores(logical=FALSE)
  if(is.na(m) || m < 1L) m <- 1L
  ans <- max(1L, as.integer(m)-1L)
  
  if(ans > 1L) {
    tm <- timecheck(dta, TS, typeTS, TSextra)
    if(length(tm) > 1L) tm <- max(tm, na.rm=TRUE)
    
    if(!is.finite(tm) || 2*tm*B < 20 || B < 2*ans) {
      ans <- 1L
      if(!SuppressMessages)
        message("maxProcessor set to 1 because parallel overhead would dominate this short computation")
    }
    else if(!SuppressMessages)
      message("Using ", ans, " cores.")
  }
  
  ans
}

#' row_named
#'
#' Row i of a matrix as a named vector; `m[i, ]` drops the name if m has one column.
#'
#' @param m a matrix
#' @param i row number
#' @return a named vector
#' @keywords internal
#' @noRd
row_named <- function(m, i) stats::setNames(m[i, ], colnames(m))

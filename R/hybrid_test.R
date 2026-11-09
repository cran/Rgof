#' Tests for the goodness-of-fit problem via twosample tests
#'
#' This function runs a number of goodness-of-fit tests in the case where data can
#' be generated under the null hypothesis but its cdf can not be calculated. A large Monte
#' Carlo data set is generated under the null hypothesis, and then the twosample tests
#' of the package R2sample are run on the data and the Monte Carlo data set.
#' p values are found via the parametric bootstrap.
#'
#' If parameters are estimated, the Monte Carlo data set is generated with the
#' parameters estimated from the data. In each bootstrap run a new data set is generated
#' with these parameters, the parameters are estimated again from this new data set,
#' and the Monte Carlo data set is generated with those estimates.
#'
#' If the cdf under the null hypothesis can be calculated, gof_test should be used
#' instead, as it generally has higher power.
#'
#' The chi square tests use the usual chi square approximation and treat the Monte
#' Carlo data as an independent sample, so their p values are only approximate if
#' parameters are estimated.
#'
#' Requires the package R2sample.
#'
#' For details on the usage of this routine consult the vignette with vignette("Rgof","Rgof")
#'
#' @param  x data set: a numeric vector for continuous data, a vector of counts for discrete data.
#' @param  vals =NA, values of discrete random variable, or NA if data is continuous.
#' @param  rnull  routine to generate data under the null hypothesis, with no arguments or
#'                one argument, the parameters.
#' @param  phat =function(x) -99, function to estimate parameters from the data, if needed.
#' @param  nMC =1, sample size of Monte Carlo data set. If nMC<=10 the sample size used
#'                 will be nMC times the sample size of x.
#' @param  TS user supplied function to find test statistics, if any. It has to have the
#'            format required by R2sample::twosample_test.
#' @param  TSextra (optional) list passed to TS, if needed.
#' @param  B =5000, number of simulation runs. If B=0 the routine returns the test statistics.
#' @param  seed =NULL (optional) seed for simulations, to aid reproducibility.
#' @param  maxProcessor number of processors to use in parallel processing.
#' @param  doMethods ="all", a vector of codes for the methods to include or all of them.
#' @param  SuppressMessages =FALSE, print informative messages?
#' @return An object of class \code{"Rgof_test"}, as returned by \code{gof_test}, with
#'   components \code{results}, \code{statistics}, \code{p.values}, \code{parameters} and
#'   \code{metadata} (which also has \code{hybrid=TRUE} and \code{nMC}, the size of the
#'   Monte Carlo data set).
#' @examples
#' # All examples are run with B=20 and maxProcessor=1 to pass CRAN checks.
#' # (they need R2sample, in a version newer than 5.0.0)
#' if(requireNamespace("R2sample", quietly=TRUE) &&
#'    utils::packageVersion("R2sample") > "5.0.0") {
#'   # Test whether data comes from a standard normal distribution.
#'   rnull=function() rnorm(100)
#'   x=rnull()
#'   hybrid_test(x, rnull=rnull, B=20, maxProcessor=1)
#'   # Test whether data comes from a normal distribution,
#'   # with mean and standard deviation estimated.
#'   rnull=function(p) rnorm(100, p[1], p[2])
#'   phat=function(x) c(mean(x), sd(x))
#'   x=rnull(c(0, 1))
#'   hybrid_test(x, rnull=rnull, phat=phat, B=20, maxProcessor=1)
#'   # Discrete data: test whether data comes from a binomial(10, p) distribution,
#'   # with p estimated
#'   vals=0:10
#'   rnull=function(p) table(c(0:10, rbinom(1000, 10, p)))-1
#'   phat=function(x) mean(rep(0:10, x))/10
#'   x=rnull(0.5)
#'   hybrid_test(x, vals, rnull, phat, B=20, maxProcessor=1)
#' }
#' @export
hybrid_test <- function(x, vals=NA, rnull, phat=function(x) -99,
                        nMC=1,
                        TS, TSextra,
                        B=5000,
                        seed=NULL,
                        maxProcessor,
                        doMethods="all",
                        SuppressMessages=FALSE) {

  check_R2sample("hybrid_test")
  Continuous <- any(is.na(vals))
  if(!Continuous) x <- as.numeric(x)
  hyb <- make_hybrid_generators(x, vals, rnull, phat, nMC, Continuous)

  # observed data and Monte Carlo data set
  if(Continuous) {
    dta <- list(x=x, y=hyb$MC(hyb$est(x)))
  } else {
    dta <- hyb$clean(list(x=x, y=hyb$MC(hyb$est(x)), vals=vals))
  }

  # Build argument list so missing optional arguments are not passed.
  args <- list(
    x = dta,
    rnull = hyb$rnullT,
    B = B,
    doMethods = doMethods,
    UseLargeSample = FALSE,
    SuppressMessages = SuppressMessages
  )
  if(!is.null(seed)) args$seed <- seed
  if(!missing(TS)) args$TS <- TS
  if(!missing(TSextra)) args$TSextra <- TSextra
  if(!missing(maxProcessor)) args$maxProcessor <- maxProcessor

  res <- do.call(R2sample::twosample_test, args)

  # same structure as gof_test, so that print and summary work the same way
  est <- hyb$est(x)
  out <- .new_Rgof_test(
    statistics=res$statistics,
    p.values=res$p.values,
    parameters=if(length(est)==1 && est[1]==-99) NULL else est,
    n=if(Continuous) length(x) else sum(x),
    Continuous=Continuous,
    B_requested=B,
    B_used=B,
    maxProcessor=if(missing(maxProcessor)) NA else maxProcessor,
    typeTS=NA,
    nbins=NA,
    rate=NA,
    minexpcount=NA,
    ChiUsePhat=NA,
    call=match.call()
  )
  out$metadata$hybrid <- TRUE
  out$metadata$nMC <- hyb$nMC
  out
}

# The hybrid tests need R2sample, in a version newer than 5.0.0: with discrete data the
# data sets of the parametric bootstrap can have a different number of values than the
# observed data, and R2sample 5.0.0 then fails with "Index out of bounds" (it uses the
# Anderson-Darling weights of the observed data for all data sets).
check_R2sample <- function(fun) {
  if(!requireNamespace("R2sample", quietly=TRUE))
    stop(fun, " requires the package R2sample", call.=FALSE)
  if(utils::packageVersion("R2sample") <= "5.0.0")
    stop(fun, " requires a version of R2sample newer than 5.0.0", call.=FALSE)
}

#' Data generators for the goodness-of-fit/twosample hybrid problem
#'
#' @param x data set
#' @param vals values of discrete random variable, or NA
#' @param rnull routine to generate data under the null hypothesis
#' @param phat routine to estimate parameters
#' @param nMC size of Monte Carlo data set, or multiple of sample size if <=10
#' @param Continuous is data continuous?
#' @keywords internal
#' @return a list of functions
#' @noRd
make_hybrid_generators <- function(x, vals, rnull, phat, nMC, Continuous) {
  nformals <- length(formals(rnull))
  gen <- function(p) {
    if(nformals==0 || (length(p)==1 && p[1]==-99)) out <- rnull()
    else out <- rnull(p)
    as.numeric(out)
  }
  est <- function(x) {
    if(nformals==0) return(-99)
    phat(x)
  }
  # size of the Monte Carlo data set
  n <- if(Continuous) length(x) else sum(x)
  if(nMC <= 10) nMC <- nMC*n
  nMC <- round(nMC)
  # Monte Carlo data set: combine data sets generated under the null hypothesis
  # until the sample size is at least nMC
  MC <- function(p) {
    y <- gen(p)
    iter <- 0
    if(Continuous) {
      while(length(y) < nMC) {
        iter <- iter+1
        if(iter>10000) stop("could not generate a Monte Carlo data set of size nMC", call.=FALSE)
        y <- c(y, gen(p))
      }
      y <- y[seq_len(nMC)]
    } else {
      while(sum(y) < nMC) {
        iter <- iter+1
        if(iter>10000) stop("could not generate a Monte Carlo data set of size nMC", call.=FALSE)
        y <- y+gen(p)
      }
    }
    y
  }
  # R2sample does not allow values with zero counts in both data sets
  clean <- function(dta) {
    I <- (dta$x+dta$y)>0
    list(x=dta$x[I], y=dta$y[I], vals=dta$vals[I])
  }
  # Parametric bootstrap, used by R2sample::twosample_test and twosample_power.
  # dta is the data set (data and Monte Carlo data). A new data set is generated
  # with the parameters estimated from the data, the parameters are estimated again
  # from the new data set, and the Monte Carlo data set is generated with these.
  rnullT <- function(dta) {
    if(Continuous) {
      xn <- gen(est(dta$x))
      return(list(x=xn, y=MC(est(xn))))
    }
    xfull <- rep(0, length(vals))
    xfull[match(dta$vals, vals)] <- dta$x
    xn <- gen(est(xfull))
    clean(list(x=xn, y=MC(est(xn)), vals=vals))
  }
  list(gen=gen, est=est, MC=MC, clean=clean, rnullT=rnullT, nMC=nMC)
}

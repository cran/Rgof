#' Power Estimation for the goodness-of-fit problem via twosample tests
#'
#' This function estimates the power of goodness-of-fit/two-sample hybrid tests
#' (see hybrid_test) by generating data under the alternative, a Monte Carlo data set
#' under the null hypothesis and then running the twosample tests of the package R2sample.
#' The null distributions of the test statistics are found via the parametric bootstrap,
#' as in hybrid_test.
#'
#' Requires the package R2sample.
#'
#' For details on the usage of this routine consult the vignette with vignette("Rgof","Rgof")
#'
#' @param  rnull  routine to generate data under the null hypothesis, with no arguments or
#'                one argument, the parameters.
#' @param  ralt  routine to generate data under the alternative hypothesis.
#' @param  param_alt values passed to ralt.
#' @param  vals =NA, values of discrete random variable, or NA if data is continuous.
#' @param  phat =function(x) -99, function to estimate parameters from the data, if needed.
#' @param  nMC =1, sample size of Monte Carlo data set. If nMC<=10 the sample size used
#'                 will be nMC times the sample size of the data.
#' @param  TS user supplied function to find test statistics, if any. It has to have the
#'            format required by R2sample::twosample_power.
#' @param  TSextra (optional) list passed to TS, if needed.
#' @param  With.p.value =FALSE, does user supplied method find its own p-values?
#' @param  alpha =0.05 type I error rate used in tests.
#' @param  B =1000  number of simulation runs.
#' @param  seed =NULL (optional) seed for simulations, to aid reproducibility.
#' @param  maxProcessor number of processors to use in parallel processing.
#' @param  doMethods ="all", a vector of codes for the methods to include or all of them.
#' @param  SuppressMessages =FALSE, print informative messages?
#' @return A matrix (or vector) of power values, with class \code{"Rgof_power_values"} so that
#'   \code{summary} gives a summary of the study (it prints and is used like the plain matrix).
#' @examples
#' # All examples are run with B=20 and maxProcessor=1 to pass CRAN checks.
#' # (they need R2sample, in a version newer than 5.0.0)
#' if(requireNamespace("R2sample", quietly=TRUE) &&
#'    utils::packageVersion("R2sample") > "5.0.0") {
#'   # Power of tests whether data comes from a standard normal distribution,
#'   # if the true distribution is normal with mean mu.
#'   rnull=function() rnorm(100)
#'   ralt=function(mu) rnorm(100, mu)
#'   hybrid_power(rnull, ralt, c(0, 0.3), B=20, maxProcessor=1)
#'   # Power of tests whether data comes from a normal distribution,
#'   # with mean and standard deviation estimated. True data comes from a t distribution.
#'   rnull=function(p) rnorm(100, p[1], p[2])
#'   ralt=function(df) rt(100, df)
#'   phat=function(x) c(mean(x), sd(x))
#'   hybrid_power(rnull, ralt, 5, phat=phat, B=20, maxProcessor=1)
#' }
#' @export
hybrid_power <- function(rnull, ralt, param_alt, vals=NA,
                         phat=function(x) -99,
                         nMC=1,
                         TS, TSextra,
                         With.p.value=FALSE,
                         alpha=0.05,
                         B=1000,
                         seed=NULL,
                         maxProcessor,
                         doMethods="all",
                         SuppressMessages=FALSE) {

  check_R2sample("hybrid_power")
  Continuous <- any(is.na(vals))
  x <- ralt(param_alt[1]) # example data set
  if(!Continuous) x <- as.numeric(x)
  hyb <- make_hybrid_generators(x, vals, rnull, phat, nMC, Continuous)

  # Function passed to R2sample::twosample_power(). For a given alternative
  # parameter a it generates the data x under the alternative and the Monte Carlo
  # data set y under the null hypothesis, with parameters estimated from x.
  f <- function(a) {
    xn <- ralt(a)
    if(Continuous) return(list(x=xn, y=hyb$MC(hyb$est(xn))))
    xn <- as.numeric(xn)
    hyb$clean(list(x=xn, y=hyb$MC(hyb$est(xn)), vals=vals))
  }

  # Build argument list so missing optional arguments are not passed.
  args <- list(
    f = f,
    param_alt,
    rnull = hyb$rnullT,
    alpha = alpha,
    With.p.value = With.p.value,
    B = B,
    SuppressMessages = SuppressMessages
  )
  if(!is.null(seed)) args$seed <- seed
  if(!missing(TS)) args$TS <- TS
  if(!missing(TSextra)) args$TSextra <- TSextra
  if(!missing(maxProcessor)) args$maxProcessor <- maxProcessor

  out <- do.call(R2sample::twosample_power, args)
  if(!(doMethods[1] %in% c("all", "default"))) {
    if(is.matrix(out)) out <- out[, doMethods, drop=FALSE]
    else out <- out[doMethods]
  }
  # plain matrix or named vector, plus the settings for summary()
  if(is.matrix(out)) out <- matrix(as.numeric(out), nrow(out), dimnames=dimnames(out))
  else out <- stats::setNames(as.numeric(out), names(out))
  n <- if(Continuous) length(x) else sum(x)
  attr(out, "metadata") <- list(
    n = n,
    data.type = if(Continuous) "continuous" else "discrete",
    hybrid = TRUE,
    nMC = if(nMC<=10) nMC*n else nMC,
    B = B,
    maxProcessor = if(missing(maxProcessor)) NA else maxProcessor,
    alpha = alpha,
    estimated = abs(phat(x)[1]+99)>0.001,
    param_alt = param_alt,
    seed = seed
  )
  class(out) <- c("Rgof_power_values", class(out))
  out
}

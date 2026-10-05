#' Tests for the univariate goodness-of-fit problem
#'
#' This function runs a number of goodness-of-fit tests using Rcpp and parallel
#' computing. The result is an object of class \code{"Rgof_test"}. The object
#' retains the traditional \code{$statistics} and \code{$p.values} components
#' and also contains a data-frame representation and metadata about the test.
#'
#' For details on the usage of this routine consult the vignette with
#' \code{vignette("Rgof", "Rgof")}.
#'
#' @param x Data set.
#' @param vals =NA, values of discrete RV, or NA if data is continuous.
#' @param pnull CDF under the null hypothesis.
#' @param rnull Routine to generate data under the null hypothesis.
#' @param w Optional function to calculate weights, returns -99 if no weights.
#' @param phat =function(x) -99, function to estimate parameters from the data,
#'   or -99 if no parameters are estimated.
#' @param TS User supplied function to find test statistics, if any.
#' @param TSextra List passed to TS, if desired, or missing.
#' @param nbins =c(50, 10), number of bins for chi-square tests.
#' @param rate =0, rate of Poisson if sample size is random, 0 if sample size is fixed.
#' @param Range =c(-Inf, Inf), limits of possible observations, if any, for
#'   chi-square tests.
#' @param B =5000, number of simulation runs. If \code{B=0}, only test
#'   statistics are calculated and p-values are returned as \code{NA}.
#' @param minexpcount =5 minimal expected bin count required.
#' @param ChiUsePhat =FALSE, if TRUE param is estimated parameter, otherwise
#'   minimum chi square method is used.
#' @param  SuppressMessages =FALSE
#' @param  seed =NULL, optional seed for reproducibility
#' @param maxProcessor Number of processors to use in parallel processing.
#' @param doMethods ="all", a vector of codes for the methods to include or all
#'   of them.
#'
#' @return An object of class \code{"Rgof_test"}. Its main components are
#'   \code{results}, a data frame with method, statistic and p-value;
#'   \code{statistics} and \code{p.values}, retained for backward
#'   compatibility; \code{parameters}, the parameter estimates when present;
#'   and \code{metadata}, information about the calculation.
#' @export
#'
#' @examples
#' pnull <- function(x) pnorm(x)
#' rnull <- function() rnorm(100)
#' x <- rnorm(100)
#' z <- gof_test(x, NA, pnull, rnull, B=500)
#' z
#' z$results
#' summary(z)
#' as.data.frame(z)
#'
gof_test <- function(x, vals=NA, pnull, rnull,
                     w=function(x) -99, phat=function(x) -99,
                     TS, TSextra, nbins=c(50, 10), rate=0,
                     Range=c(-Inf, Inf), B=5000, minexpcount=5.0,
                     ChiUsePhat = FALSE, SuppressMessages=FALSE,
                     seed=NULL, maxProcessor, doMethods="all") {

  call <- match.call()
  B_requested <- B

  if(any(is.na(vals))) {
    Continuous <- TRUE
    dta <- list(x=x)
    check.functions(pnull, rnull, phat, x=x)
  }
  else {
    Continuous <- FALSE
    dta <- list(x=x, vals=vals)
    check.functions(pnull, rnull, phat, vals, x)
  }

  # Estimate once for metadata and for decisions that depend on whether
  # parameters are estimated.
  phat_data <- phat(x)
  has_estimated_parameters <- !(length(phat_data) == 0L ||
                                abs(phat_data[1] + 99) < 0.001)
  parameters <- if(has_estimated_parameters) phat_data else numeric(0)

  TSextra <- makeTSextra(TSextra, pnull, phat, w, Continuous)

  WithWeights <- TRUE
  if(length(formals(w)) == 1L) {
    if(w(x[1]) == -99) WithWeights <- FALSE
  }

  if(test_methods(doMethods, Continuous, WithWeights)) return(NULL)

  # Adjust number of bins to account for parameter estimation.
  if(has_estimated_parameters)
    nbins <- nbins + length(phat_data)

  tmp <- maketypeTS(TS, Continuous, WithWeights)
  typeTS <- tmp$typeTS
  TS <- tmp$TS
  useSingleProcessor <- tmp$useSingleProcessor
  if(!is.null(seed)) {
    if(length(seed)!=1L || is.na(seed) || !is.finite(seed))
      stop("seed must be NULL or a single finite number", call.=FALSE)
    
    set.seed(seed)
  }
  TS_data <- calcTS(dta, TS, typeTS, TSextra)
  if(is.null(names(TS_data))) {
    message("result of TS has to be a named vector")
    return(NULL)
  }

  # B=0: retain the structured return type, but no simulated p-values exist.
  if(B == 0) {
    statistics <- TS_data
    p.values <- rep(NA_real_, length(statistics))
    names(p.values) <- names(statistics)

    return(.new_Rgof_test(
      statistics=statistics,
      p.values=p.values,
      parameters=parameters,
      n=if(Continuous) length(x) else sum(x),
      Continuous=Continuous,
      B_requested=B_requested,
      B_used=0L,
      maxProcessor=NA_integer_,
      typeTS=typeTS,
      nbins=nbins,
      rate=rate,
      minexpcount=minexpcount,
      ChiUsePhat=ChiUsePhat,
      call=call
    ))
  }

  # Choose number of processors after dta, TS and typeTS are known.
  maxProcessor <- makemaxProcessor(maxProcessor,
      dta=dta, TS=TS, typeTS=typeTS, TSextra=TSextra, B=B,
      SuppressMessages=SuppressMessages,
      useSingleProcessor=useSingleProcessor
    )
  # Assure that B is a multiple of maxProcessor.
  B=ceiling(B/maxProcessor)*maxProcessor

  if(length(x) > 10000 && maxProcessor == 1)
    message("Consider using parallel processing with maxProcessor= (your number of cores)")

  if(maxProcessor == 1) {
    outTS <- gof_test_C(dta, rnull, TS, typeTS, TSextra, B)
  }
  else {
    cl <- parallel::makeCluster(maxProcessor)
    on.exit(parallel::stopCluster(cl), add=TRUE)
    if(!is.null(seed))
      parallel::clusterSetRNGStream(cl, iseed=seed)
    z <- parallel::clusterCall(
      cl, gof_test_C,
      dta, rnull, TS, typeTS, TSextra,
      B=round(B/maxProcessor)
    )

    # Average results across cores.
    tmp <- 0*z[[1]]
    for(i in seq_len(maxProcessor)) tmp <- tmp + z[[i]]
    outTS <- tmp/maxProcessor
  }

  # Do chi-square tests.
  if(typeTS == 1) {
    if(is.infinite(Range[1])) Range[1] <- -99999
    if(is.infinite(Range[2])) Range[2] <- 99999
    outchi <- t(chi_test_cont(
      x, pnull, w, phat,
      ifelse(TSextra$Noqnull, NA, TSextra$qnull),
      nbins, rate, Range, minexpcount, ChiUsePhat
    )[, c(1, 2)])
  }

  if(typeTS == 5) {
    outchi <- t(chi_test_disc(
      x, pnull, phat,
      nbins, rate, minexpcount, ChiUsePhat
    )[, 1:2])
  }

  if(typeTS == 1 || typeTS == 5) {
    out <- list(
      statistics=c(outTS[1, ], outchi[1, ]),
      p.values=c(outTS[2, ], outchi[2, ])
    )
  }
  else {
    out <- list(statistics=outTS[1, ], p.values=outTS[2, ])
  }

  if(doMethods[1] != "all") {
    out[[1]] <- out[[1]][doMethods]
    out[[2]] <- out[[2]][doMethods]
  }

  if(TSextra$Noqnull) {
    out[[1]] <- out[[1]][names(out[[1]]) != "Wassp1"]
    out[[2]] <- out[[2]][names(out[[2]]) != "Wassp1"]
  }

  # Keep the existing presentation rounding.
  out <- signif_digits(out)

  .new_Rgof_test(
    statistics=out$statistics,
    p.values=out$p.values,
    parameters=parameters,
    n=if(Continuous) length(x) else sum(x),
    Continuous=Continuous,
    B_requested=B_requested,
    B_used=B,
    maxProcessor=maxProcessor,
    typeTS=typeTS,
    nbins=nbins,
    rate=rate,
    minexpcount=minexpcount,
    ChiUsePhat=ChiUsePhat,
    call=call
  )
}


# Internal constructor for structured gof_test output.
.new_Rgof_test <- function(statistics, p.values, parameters, n, Continuous,
                           B_requested, B_used, maxProcessor, typeTS, nbins,
                           rate, minexpcount, ChiUsePhat, call) {

  methods <- names(statistics)
  if(is.null(methods)) methods <- rep("", length(statistics))

  if(length(p.values) != length(statistics))
    stop("statistics and p.values must have the same length", call.=FALSE)

  results <- data.frame(
    method=methods,
    statistic=unname(statistics),
    p.value=unname(p.values),
    row.names=NULL,
    check.names=FALSE
  )

  structure(
    list(
      results=results,
      statistics=statistics,
      p.values=p.values,
      parameters=parameters,
      metadata=list(
        n=n,
        data.type=if(Continuous) "continuous" else "discrete",
        B.requested=B_requested,
        B.used=B_used,
        maxProcessor=maxProcessor,
        typeTS=typeTS,
        nbins=nbins,
        rate=rate,
        minexpcount=minexpcount,
        ChiUsePhat=ChiUsePhat
      ),
      call=call
    ),
    class="Rgof_test"
  )
}


#' Print an Rgof goodness-of-fit test result
#'
#' @param x An object of class \code{"Rgof_test"}.
#' @param ... Further arguments passed to \code{print.data.frame}.
#' @return \code{x}, invisibly.
#' @export
print.Rgof_test <- function(x, ...) {
  cat("\nRgof goodness-of-fit tests\n")
  cat("Data:", x$metadata$data.type,
      " | n =", x$metadata$n,
      " | B =", x$metadata$B.used, "\n\n")

  print(x$results, row.names=FALSE, ...)

  if(length(x$parameters)) {
    cat("\nEstimated parameters:\n")
    print(x$parameters)
  }

  invisible(x)
}


#' Summarize an Rgof goodness-of-fit test result
#'
#' @param object An object of class \code{"Rgof_test"}.
#' @param ... Currently unused.
#' @return An object of class \code{"summary.Rgof_test"}.
#' @export
summary.Rgof_test <- function(object, ...) {
  out <- list(
    call=object$call,
    results=object$results,
    parameters=object$parameters,
    metadata=object$metadata
  )
  class(out) <- "summary.Rgof_test"
  out
}


#' @export
print.summary.Rgof_test <- function(x, ...) {
  cat("\nRgof goodness-of-fit test summary\n\n")
  cat("Data type: ", x$metadata$data.type, "\n", sep="")
  cat("Sample size: ", x$metadata$n, "\n", sep="")
  cat("Simulation runs requested: ", x$metadata$B.requested, "\n", sep="")
  cat("Simulation runs used: ", x$metadata$B.used, "\n", sep="")

  if(!is.na(x$metadata$maxProcessor))
    cat("Processors used: ", x$metadata$maxProcessor, "\n", sep="")

  if(length(x$parameters)) {
    cat("\nEstimated parameters:\n")
    print(x$parameters)
  }

  cat("\nTest results:\n")
  print(x$results, row.names=FALSE)
  invisible(x)
}


#' Convert an Rgof test result to a data frame
#'
#' @param x An object of class \code{"Rgof_test"}.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Currently unused.
#' @return A data frame with columns \code{method}, \code{statistic}, and
#'   \code{p.value}.
#' @export
as.data.frame.Rgof_test <- function(x, row.names=NULL, optional=FALSE, ...) {
  x$results
}

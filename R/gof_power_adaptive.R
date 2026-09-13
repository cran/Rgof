#' Adaptive power estimation for goodness-of-fit tests
#'
#' Estimates power in simulation batches and stops when the Monte Carlo
#' standard error of every reported power estimate is no larger than
#' \code{target.mcse}, or when \code{B.max} is reached.
#'
#' Unlike \code{gof_power()}, this routine performs the adaptive simulation
#' directly.  For statistic-based tests it estimates the null critical values
#' once using \code{B.null} null simulations and holds those critical values
#' fixed while alternative simulations are added in batches.
#'
#' @param pnull Function to calculate the cdf under the null hypothesis.
#' @param vals =NA values of a discrete random variable, or NA for continuous data.
#' @param rnull Function to generate data under the null hypothesis.
#' @param ralt Function to generate data under the alternative hypothesis.
#' @param param_alt Vector of parameter values under the alternative.
#' @param w Optional weight function; returns -99 if no weights are used.
#' @param phat Function to estimate parameters, or function(x) -99.
#' @param TS Optional user-supplied test statistic or p-value routine.
#' @param TSextra Optional list supplied to TS.
#' @param With.p.value =FALSE; TRUE if a user-supplied TS returns p-values.
#' @param alpha Significance level.
#' @param Range Limits of possible continuous observations.
#' @param nbins Number of bins for chi-square tests.
#' @param rate Poisson rate if sample size is random; 0 for fixed sample size.
#' @param maxProcessor Maximum number of processors used for the statistic-based
#'   simulation.  P-value and chi-square simulation remains sequential, as in
#'   the corresponding existing package routines.
#' @param minexpcount Minimum expected bin count for chi-square tests.
#' @param ChiUsePhat If TRUE, use estimated parameters in chi-square tests.
#' @param SuppressMessages Suppress adaptive progress messages when TRUE.
#'   A run-time estimate is still displayed when the estimated maximum run
#'   time exceeds 30 seconds.
#' @param target.mcse Target Monte Carlo standard error.
#' @param B.min Minimum number of alternative simulations per parameter.
#' @param B.max Maximum number of alternative simulations per parameter.
#' @param B.batch Number of additional alternative simulations per batch.
#' @param B.null Number of null simulations used once to estimate the fixed
#'   critical values for statistic-based tests. This is independent of
#'   \code{B.min} and \code{B.max}, which control alternative simulations.
#' @param conf.level Confidence level for Wilson confidence intervals.
#' @param list.with.everything Optional case-study list accepted by the package.
#'
#' @return An object of class \code{Rgof_power_adaptive} and \code{Rgof_power}.
#'   It contains power estimates, Monte Carlo standard errors, Wilson confidence
#'   intervals, total alternative simulation count, null simulation count,
#'   convergence information, and (when applicable) the fixed critical values.
#'
#' @details
#' Statistic-based tests use a two-stage procedure. First, \code{B.null} null
#' data sets are simulated and one set of critical values is estimated. Those
#' critical values are then held fixed. The alternative simulations generated
#' incidentally by that \code{gof_power_C()} call are discarded. Alternative
#' power simulation then starts separately at \code{B.min} and proceeds in
#' batches, accumulating exact rejection counts up to \code{B.max}.
#'
#' User tests with \code{With.p.value=TRUE} do not require null critical-value
#' simulation; rejection is determined directly from p-value < alpha.
#'
#' The package's chi-square tests already perform their own bin construction and
#' diagnostics. Their rejection counts are accumulated directly from
#' \code{chi_test_cont()} or \code{chi_test_disc()}. For discrete data the bin
#' definitions are selected once for each alternative parameter, matching the
#' design of \code{chi_power_disc()}.
#'
#' @export
gof_power_adaptive <- function(
    pnull, vals=NA, rnull, ralt, param_alt,
    w=function(x) -99, phat=function(x) -99, TS, TSextra,
    With.p.value=FALSE,
    alpha=0.05, Range=c(-Inf, Inf), nbins=c(50, 10),
    rate=0, maxProcessor, minexpcount=5.0, ChiUsePhat=TRUE,
    SuppressMessages=FALSE,
    target.mcse=0.01,
    B.min=500,
    B.max=10000,
    B.batch=250,
    B.null=1000,
    conf.level=0.95,
    list.with.everything) {

  fff <- nortest::lillie.test # avoid CRAN namespace issues; intentionally unused

  ## ---- validate adaptive controls ---------------------------------------
  scalar_integer <- function(x, nm, minimum=1L) {
    if(!is.numeric(x) || length(x)!=1L || is.na(x) ||
       !is.finite(x) || x!=floor(x) || x<minimum)
      stop(paste0(nm, " must be an integer >= ", minimum), call.=FALSE)
    as.integer(x)
  }

  if(!is.numeric(target.mcse) || length(target.mcse)!=1L ||
     is.na(target.mcse) || !is.finite(target.mcse) ||
     target.mcse<=0 || target.mcse>=0.5)
    stop("target.mcse must be a single number between 0 and 0.5",
         call.=FALSE)

  B.min <- scalar_integer(B.min, "B.min")
  B.max <- scalar_integer(B.max, "B.max")
  B.batch <- scalar_integer(B.batch, "B.batch")
  B.null <- scalar_integer(B.null, "B.null")

  if(B.min>B.max)
    stop("B.min must not exceed B.max", call.=FALSE)

  if(!is.numeric(conf.level) || length(conf.level)!=1L ||
     is.na(conf.level) || !is.finite(conf.level) ||
     conf.level<=0 || conf.level>=1)
    stop("conf.level must be a single number between 0 and 1",
         call.=FALSE)

  if(!is.numeric(alpha) || length(alpha)!=1L || is.na(alpha) ||
     alpha<=0 || alpha>=1)
    stop("alpha must be a single number between 0 and 1",
         call.=FALSE)

  ## ---- case-study input --------------------------------------------------
  if(!missing(list.with.everything)) {
    pnull <- list.with.everything$pnull
    vals <- list.with.everything$vals
    rnull <- list.with.everything$rnull

    if(missing(ralt))
      ralt <- list.with.everything$ralt

    phat <- list.with.everything$phat

    if(!is.null(list.with.everything$Range))
      Range <- list.with.everything$Range

    case_extra <- list.with.everything$TSextra
    if(is.null(case_extra))
      case_extra <- list()

    if(missing(TSextra))
      TSextra <- case_extra
    else
      TSextra <- utils::modifyList(case_extra, TSextra)
  }

  NewTest <- !missing(TS)

  ## ralt() with no formal arguments is converted to the one-argument
  ## convention expected throughout the power code.
  if(length(formals(ralt))==0L) {
    ralt0 <- ralt
    ralt <- function(n) ralt0()
    param_alt <- 0
  }

  ## Example data set, as in gof_power(), used for validation and dispatch.
  x <- ralt(param_alt[1])
  Continuous <- any(is.na(vals))

  if(Continuous) {
    dta <- list(x=x)
    check.functions(pnull, rnull, phat, x=x)
  } else {
    dta <- list(x=x, vals=vals)
    check.functions(pnull, rnull, phat, vals, x)
  }

  TSextra <- makeTSextra(TSextra, pnull, phat, w, Continuous)

  WithWeights <- TRUE
  if(length(formals(w))==1L) {
    if(w(x[1])==-99)
      WithWeights <- FALSE
  }

  ## Account for estimated parameters in the chi-square bin choices.
  if(abs(phat(x)[1]+99)>0.001)
    nbins <- nbins+length(phat(x))

  tmp <- maketypeTS(TS, Continuous, WithWeights)
  typeTS <- tmp$typeTS
  TS <- tmp$TS

  TS_data <- calcTS(dta, TS, typeTS, TSextra)
  if(is.null(names(TS_data)))
    stop("result of TS has to be a named vector", call.=FALSE)

  ## ---- processor setup ---------------------------------------------------
  cores <- parallel::detectCores(logical=FALSE)
  if(is.na(cores))
    cores <- 1L

  available <- max(1L, cores-1L)

  if(missing(maxProcessor))
    maxProcessor <- available

  maxProcessor <- scalar_integer(maxProcessor, "maxProcessor")
  maxProcessor <- min(maxProcessor, available)

  ## Existing p-value implementation is sequential.
  if(With.p.value)
    maxProcessor <- 1L

  ## Custom C/C++ TS cannot safely use the package's PSOCK strategy.
  if(NewTest) {
    ts_txt <- tryCatch(deparse(TS), error=function(e) character())
    if(any(grepl("\\.Call", ts_txt, fixed=FALSE))) {
      if(maxProcessor>1L && !SuppressMessages)
        message("Parallel programming is not possible for this custom C/C++ TS; using one processor")
      maxProcessor <- 1L
    }
  }

  ## For short runs, follow the same time-based policy as gof_power().
  if(maxProcessor>1L) {
    tm <- timecheck(dta, TS, typeTS, TSextra)
    if(tm*length(param_alt)*B.min<20) {
      maxProcessor <- 1L
      if(!SuppressMessages)
        message("maxProcessor set to 1 for faster computation")
    } else if(!SuppressMessages) {
      message(paste("Using", maxProcessor, "cores.."))
    }
  }

  ## ---- early run-time estimate ------------------------------------------
  ## Estimate total run time from the same one-data-set timing probe already
  ## used by the package.  Display it only when the estimated run time exceeds
  ## 30 seconds.  The estimate is deliberately based on B.max, so it is an
  ## upper-bound-style estimate; adaptive stopping may finish substantially
  ## earlier.
  if(!With.p.value) {
    tm.estimate <- timecheck(dta, TS, typeTS, TSextra)
    estimated.seconds <- tm.estimate *
      (B.null + length(param_alt)*B.max) / maxProcessor

    if(typeTS==1L || (typeTS==5L && !NewTest))
      estimated.seconds <- 2*estimated.seconds

    if(is.finite(estimated.seconds) && estimated.seconds>30) {
      if(estimated.seconds<60) {
        runtime.text <- sprintf("%.0f seconds", estimated.seconds)
      } else if(estimated.seconds<3600) {
        runtime.text <- sprintf("%.1f minutes", estimated.seconds/60)
      } else {
        runtime.text <- sprintf("%.1f hours", estimated.seconds/3600)
      }

      message(
        "Estimated maximum run time: ", runtime.text,
        ". Adaptive stopping may finish earlier."
      )
    }
  }

  ## ---- helpers -----------------------------------------------------------

  ## Run gof_power_C directly for one exact-size statistic batch.
  ## Unlike gof_power(), this helper distributes an arbitrary B exactly across
  ## workers, so adaptive batches do not have to be multiples of maxProcessor.
  run_ts_batch <- function(Bnow) {

    if(maxProcessor==1L || Bnow==1L)
      return(gof_power_C(rnull, vals, ralt, param_alt,
                         TS, typeTS, TSextra, Bnow))

    nw <- min(maxProcessor, Bnow)
    reps <- rep(Bnow %/% nw, nw)
    if(Bnow %% nw)
      reps[seq_len(Bnow %% nw)] <- reps[seq_len(Bnow %% nw)]+1L

    cl <- parallel::makeCluster(nw)
    on.exit(parallel::stopCluster(cl), add=TRUE)

    z <- parallel::clusterMap(
      cl,
      gof_power_C,
      B=reps,
      MoreArgs=list(
        rnull=rnull,
        vals=vals,
        ralt=ralt,
        param_alt=param_alt,
        TS=TS,
        typeTS=typeTS,
        TSextra=TSextra
      ),
      SIMPLIFY=FALSE,
      USE.NAMES=FALSE
    )

    Data <- do.call(rbind, lapply(z, `[[`, "Data"))
    Sim <- do.call(rbind, lapply(z, `[[`, "Sim"))
    list(Data=Data, Sim=Sim)
  }

  ## Convert a Sim matrix into exact rejection counts using fixed critical
  ## values. Rows of Sim are grouped by alternative parameter.
  count_ts_rejections <- function(Sim, critical.values) {
    out <- matrix(0L, length(param_alt), length(critical.values),
                  dimnames=list(as.character(param_alt), names(TS_data)))

    for(i in seq_along(param_alt)) {
      ss <- Sim[Sim[,1]==param_alt[i], -1, drop=FALSE]

      for(j in seq_along(critical.values))
        out[i,j] <- sum(ss[,j]>critical.values[j])
    }
    out
  }

  ## Direct simulation for a custom routine that returns p-values.
  count_pvalue_rejections <- function(Bnow) {
    out <- matrix(0L, length(param_alt), length(TS_data),
                  dimnames=list(as.character(param_alt), names(TS_data)))

    dd <- dta

    for(i in seq_along(param_alt)) {
      for(b in seq_len(Bnow)) {
        dd$x <- ralt(param_alt[i])
        pv <- calcTS(dd, TS, typeTS, TSextra)
        out[i,] <- out[i,] + as.integer(pv<alpha)
      }
    }
    out
  }

  ## Continuous chi-square methods. chi_test_cont() performs the package's
  ## normal bin construction/diagnostics for each generated data set.
  count_chi_cont_rejections <- function(Bnow) {
    method_names <- c("ES-l-P", "ES-s-P", "EP-l-P", "EP-s-P",
                      "ES-l-L", "ES-s-L", "EP-l-L", "EP-s-L")

    out <- matrix(0L, length(param_alt), length(method_names),
                  dimnames=list(as.character(param_alt), method_names))

    qfun <- if(TSextra$Noqnull) NA else TSextra$qnull

    RR <- Range
    if(is.infinite(RR[1])) RR[1] <- -99999
    if(is.infinite(RR[2])) RR[2] <- 99999

    for(i in seq_along(param_alt)) {
      for(b in seq_len(Bnow)) {
        xx <- ralt(param_alt[i])

        pv <- chi_test_cont(
          xx, pnull, w, phat, qfun,
          nbins, rate, RR, minexpcount, ChiUsePhat
        )[,2]

        out[i,] <- out[i,] + as.integer(pv<alpha)
      }
    }
    out
  }

  ## For discrete data the existing chi_power_disc() chooses its bin groups
  ## once per alternative parameter. Do the same here, then retain those bins
  ## for all adaptive batches.
  disc_bins <- NULL
  if(typeTS==5L && !NewTest) {
    disc_bins <- lapply(seq_along(param_alt), function(i) {
      xx <- ralt(param_alt[i])
      make_bins_disc(xx, pnull, phat,
                     nbins=nbins, minexpcount=minexpcount)
    })
  }

  count_chi_disc_rejections <- function(Bnow) {
    ## gof_power() currently keeps the first two discrete chi-square powers.
    method_names <- c("l-P", "s-P")

    out <- matrix(0L, length(param_alt), length(method_names),
                  dimnames=list(as.character(param_alt), method_names))

    for(i in seq_along(param_alt)) {
      for(b in seq_len(Bnow)) {
        xx <- ralt(param_alt[i])

        pv <- chi_test_disc(
          xx, pnull, phat,
          rate=rate,
          minexpcount=minexpcount,
          ChiUsePhat=ChiUsePhat,
          allbins=disc_bins[[i]]
        )[,2]

        out[i,] <- out[i,] + as.integer(pv[1:2]<alpha)
      }
    }
    out
  }

  combine_counts <- function(main, chi) {
    if(is.null(chi))
      main
    else
      cbind(main, chi)
  }

  ## ---- stage 1: fixed null critical values -------------------------------
  critical.values <- NULL
  main.count <- NULL
  B.used <- 0L

  if(With.p.value) {

    first.B <- min(B.min, B.max)
    main.count <- count_pvalue_rejections(first.B)
    B.used <- first.B

  } else {

    ## Stage 1: estimate null critical values once.
    ##
    ## gof_power_C() also returns alternative simulations, but those are
    ## intentionally discarded here. B.null controls only the null critical-
    ## value stage; B.min/B.max control the alternative power stage.
    null.run <- run_ts_batch(B.null)

    critical.values <- apply(null.run$Data, 2, stats::quantile,
                             prob=1-alpha, na.rm=TRUE)
    names(critical.values) <- names(TS_data)

    ## Stage 2: begin alternative simulation independently of B.null.
    ## run_ts_batch() necessarily generates another null matrix internally,
    ## because that is the interface of gof_power_C(); those null rows are
    ## ignored after the fixed critical values have been obtained.
    first.B <- min(B.min, B.max)
    alt.run <- run_ts_batch(first.B)

    main.count <- count_ts_rejections(alt.run$Sim, critical.values)
    B.used <- first.B
  }

  ## Chi-square rejection counts use separate alternative samples, matching
  ## the existing gof_power()/chi_power_* design.
  chi.count <- NULL

  if(typeTS==1L) {
    chi.count <- count_chi_cont_rejections(B.used)
  }

  if(typeTS==5L && !NewTest) {
    chi.count <- count_chi_disc_rejections(B.used)
  }

  rejection.count <- combine_counts(main.count, chi.count)

  ## ---- adaptive alternative simulation ----------------------------------
  converged <- FALSE

  repeat {

    power <- rejection.count/B.used
    mc.se <- sqrt(power*(1-power)/B.used)

    if(B.used>=B.min &&
       all(is.finite(mc.se)) &&
       all(mc.se<=target.mcse)) {
      converged <- TRUE
      break
    }

    if(B.used>=B.max)
      break

    this.B <- min(B.batch, B.max-B.used)

    if(With.p.value) {
      main.new <- count_pvalue_rejections(this.B)
    } else {
      nxt <- run_ts_batch(this.B)
      main.new <- count_ts_rejections(nxt$Sim, critical.values)
    }

    chi.new <- NULL

    if(typeTS==1L)
      chi.new <- count_chi_cont_rejections(this.B)

    if(typeTS==5L && !NewTest)
      chi.new <- count_chi_disc_rejections(this.B)

    main.count <- main.count+main.new

    if(!is.null(chi.new))
      chi.count <- chi.count+chi.new

    rejection.count <- combine_counts(main.count, chi.count)
    B.used <- B.used+this.B

    if(!SuppressMessages) {
      current.max <- max(
        sqrt((rejection.count/B.used)*
               (1-rejection.count/B.used)/B.used),
        na.rm=TRUE
      )
      message(sprintf("B = %d, maximum MCSE = %.5f",
                      B.used, current.max))
    }
  }

  power <- rejection.count/B.used
  mc.se <- sqrt(power*(1-power)/B.used)

  ## ---- Wilson intervals --------------------------------------------------
  z <- stats::qnorm(1-(1-conf.level)/2)
  den <- 1+z^2/B.used
  center <- (power+z^2/(2*B.used))/den
  half <- z*sqrt(power*(1-power)/B.used +
                   z^2/(4*B.used^2))/den

  ## Matrix expression first: preserve dimensions and dimnames.
  lower <- pmax(center-half, 0)
  upper <- pmin(center+half, 1)

  ## Match gof_power() convention: simplify one-alternative matrices to
  ## named vectors, but explicitly retain method names.
  if(is.matrix(power) && nrow(power)==1L) {
    method_names <- colnames(power)

    power <- as.numeric(power[1,])
    mc.se <- as.numeric(mc.se[1,])
    lower <- as.numeric(lower[1,])
    upper <- as.numeric(upper[1,])
    rejection.count <- as.numeric(rejection.count[1,])

    names(power) <- method_names
    names(mc.se) <- method_names
    names(lower) <- method_names
    names(upper) <- method_names
    names(rejection.count) <- method_names
  }

  if(!converged && !SuppressMessages)
    message("B.max reached before the target Monte Carlo SE was achieved")

  out <- list(
    power=power,
    mc.se=mc.se,
    lower=lower,
    upper=upper,
    rejections=rejection.count,
    B=B.used,
    B.null=if(With.p.value) 0L else B.null,
    critical.values=critical.values,
    target.mcse=target.mcse,
    converged=converged,
    B.min=B.min,
    B.max=B.max,
    B.batch=B.batch,
    conf.level=conf.level,
    interval="Wilson",
    alpha=alpha,
    param_alt=param_alt,
    call=match.call()
  )

  class(out) <- c("Rgof_power_adaptive", "Rgof_power")
  out
}


#' Print an adaptive Rgof power object
#'
#' @param x Object returned by \code{gof_power_adaptive()}.
#' @param ... Unused.
#' @export
print.Rgof_power_adaptive <- function(x, ...) {

  cat("Adaptive Monte Carlo power estimation\n")
  cat("Alternative simulations:", x$B,
      "| null simulations:", x$B.null,
      "| target MCSE:", x$target.mcse,
      "| converged:", if(x$converged) "yes" else "no",
      "\n")
  cat(sprintf("%.1f%% Wilson confidence intervals\n",
              100*x$conf.level))

  if(is.matrix(x$power)) {
    for(i in seq_len(nrow(x$power))) {
      if(!is.null(rownames(x$power)))
        cat("\nAlternative parameter:", rownames(x$power)[i], "\n")

      tab <- data.frame(
        power=x$power[i,],
        mc.se=x$mc.se[i,],
        lower=x$lower[i,],
        upper=x$upper[i,],
        rejections=x$rejections[i,],
        check.names=FALSE
      )

      if(!is.null(colnames(x$power)))
        rownames(tab) <- colnames(x$power)

      print(tab)
    }
  } else {
    tab <- data.frame(
      power=x$power,
      mc.se=x$mc.se,
      lower=x$lower,
      upper=x$upper,
      rejections=x$rejections,
      check.names=FALSE
    )

    if(!is.null(names(x$power)))
      rownames(tab) <- names(x$power)

    print(tab)
  }

  invisible(x)
}


#' Convert adaptive Rgof power output to a data frame
#'
#' @param x Object returned by \code{gof_power_adaptive()}.
#' @param ... Unused.
#' @export
as.data.frame.Rgof_power_adaptive <- function(x, ...) {

  if(is.matrix(x$power)) {
    nr <- nrow(x$power)
    nc <- ncol(x$power)

    alt <- rownames(x$power)
    if(is.null(alt))
      alt <- as.character(seq_len(nr))

    method <- colnames(x$power)
    if(is.null(method))
      method <- as.character(seq_len(nc))

    return(data.frame(
      param_alt=rep(alt, times=nc),
      method=rep(method, each=nr),
      power=c(x$power),
      mc.se=c(x$mc.se),
      lower=c(x$lower),
      upper=c(x$upper),
      rejections=c(x$rejections),
      row.names=NULL,
      check.names=FALSE
    ))
  }

  method <- names(x$power)
  if(is.null(method))
    method <- as.character(seq_along(x$power))

  data.frame(
    param_alt=rep(as.character(x$param_alt[1]), length(x$power)),
    method=method,
    power=unname(x$power),
    mc.se=unname(x$mc.se),
    lower=unname(x$lower),
    upper=unname(x$upper),
    rejections=unname(x$rejections),
    row.names=NULL,
    check.names=FALSE
  )
}

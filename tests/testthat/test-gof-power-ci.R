local_power_ci_setup <- function(n = 40) {
  list(
    pnull = function(x) pnorm(x),
    rnull = function() rnorm(n),
    ralt = function(mu) rnorm(n, mu),
    TSextra = list(qnull = function(x) qnorm(x))
  )
}

test_that("gof_power default return remains numeric", {
  set.seed(321)
  z <- local_power_ci_setup()

  out <- gof_power(
    z$pnull, NA, z$rnull, z$ralt,
    param_alt = c(0, 0.5),
    TSextra = z$TSextra,
    B = 12, maxProcessor = 1
  )

  expect_true(is.matrix(out))
  expect_false(inherits(out, "Rgof_power"))
  expect_true(all(out >= 0 & out <= 1))
})

test_that("gof_power returns Wilson confidence intervals when requested", {
  set.seed(322)
  z <- local_power_ci_setup()

  myTS <- function(x, pnull, param) {
    out <- abs(mean(x))
    names(out) <- "MeanAbs"
    out
  }

  out <- gof_power(
    z$pnull, NA, z$rnull, z$ralt,
    param_alt = c(0, 0.5),
    TS = myTS,
    B = 20, maxProcessor = 1,
    CI = TRUE
  )

  expect_s3_class(out, "Rgof_power")
  expect_named(
    out,
    c("power", "mc.se", "lower", "upper", "B",
      "conf.level", "interval", "alpha", "param_alt", "seed", "metadata")
  )
  expect_type(out$metadata, "list")
  expect_equal(out$metadata$B, 20)

  expect_true(is.matrix(out$power))
  expect_equal(dim(out$power), c(2L, 1L))
  expect_equal(dim(out$lower), dim(out$power))
  expect_equal(dim(out$upper), dim(out$power))
  expect_equal(dim(out$mc.se), dim(out$power))

  expect_equal(out$B, 20)
  expect_equal(out$conf.level, 0.95)
  expect_identical(out$interval, "Wilson")

  expect_true(all(out$lower >= 0 & out$lower <= 1))
  expect_true(all(out$upper >= 0 & out$upper <= 1))
  expect_true(all(out$lower <= out$power))
  expect_true(all(out$power <= out$upper))
  expect_true(all(out$mc.se >= 0))
})

test_that("gof_power supports custom confidence levels", {
  set.seed(323)
  z <- local_power_ci_setup()

  myTS <- function(x, pnull, param) {
    out <- abs(mean(x))
    names(out) <- "MeanAbs"
    out
  }

  out90 <- gof_power(
    z$pnull, NA, z$rnull, z$ralt,
    param_alt = c(0, 0.5), TS = myTS,
    B = 20, maxProcessor = 1,
    CI = TRUE, conf.level = 0.90
  )

  expect_equal(out90$conf.level, 0.90)
  expect_true(all(out90$lower <= out90$power))
  expect_true(all(out90$power <= out90$upper))
})

test_that("gof_power CI output handles one alternative", {
  set.seed(324)
  z <- local_power_ci_setup()

  myTS <- function(x, pnull, param) {
    out <- abs(mean(x))
    names(out) <- "MeanAbs"
    out
  }

  out <- gof_power(
    z$pnull, NA, z$rnull, z$ralt,
    param_alt = 0.5, TS = myTS,
    B = 20, maxProcessor = 1,
    CI = TRUE
  )

  expect_s3_class(out, "Rgof_power")
  expect_true(is.numeric(out$power))
  expect_false(is.matrix(out$power))
  expect_named(out$power, "MeanAbs")
  expect_named(out$lower, "MeanAbs")
  expect_named(out$upper, "MeanAbs")
})

test_that("as.data.frame works for Rgof_power confidence intervals", {
  set.seed(325)
  z <- local_power_ci_setup()

  myTS <- function(x, pnull, param) {
    out <- abs(mean(x))
    names(out) <- "MeanAbs"
    out
  }

  out <- gof_power(
    z$pnull, NA, z$rnull, z$ralt,
    param_alt = c(0, 0.5), TS = myTS,
    B = 20, maxProcessor = 1,
    CI = TRUE
  )

  d <- as.data.frame(out)

  expect_s3_class(d, "data.frame")
  expect_named(d, c("param_alt", "method", "power", "mc.se", "lower", "upper"))
  expect_equal(nrow(d), length(out$power))
  expect_true(all(d$lower <= d$power))
  expect_true(all(d$power <= d$upper))
})

test_that("gof_power validates CI arguments", {
  z <- local_power_ci_setup()

  expect_error(
    gof_power(
      z$pnull, NA, z$rnull, z$ralt,
      param_alt = 0, TSextra = z$TSextra,
      B = 8, maxProcessor = 1,
      CI = NA
    ),
    "CI must be TRUE or FALSE"
  )

  expect_error(
    gof_power(
      z$pnull, NA, z$rnull, z$ralt,
      param_alt = 0, TSextra = z$TSextra,
      B = 8, maxProcessor = 1,
      CI = TRUE, conf.level = 1
    ),
    "conf.level must be a single number between 0 and 1"
  )
})

test_that("summary of gof_power gives a study summary", {
  set.seed(321)
  z <- local_power_ci_setup()

  out <- gof_power(
    z$pnull, NA, z$rnull, z$ralt,
    param_alt = c(0, 0.5),
    TSextra = z$TSextra,
    B = 12, maxProcessor = 1
  )
  s <- summary(out)
  expect_s3_class(s, "summary.Rgof_power")
  expect_equal(unname(s$power), unname(unclass(out)[, , drop = FALSE]),
               ignore_attr = TRUE)
  expect_equal(s$metadata$data.type, "continuous")
  expect_equal(s$metadata$n, 40)
  expect_output(print(s), "Rgof power study summary")
  expect_output(print(out), "KS")
  # subsetting gives plain numbers
  expect_false(inherits(out[1, ], "Rgof_power_values"))

  ci <- gof_power(
    z$pnull, NA, z$rnull, z$ralt,
    param_alt = 0.5,
    TSextra = z$TSextra,
    B = 12, maxProcessor = 1, CI = TRUE
  )
  expect_output(print(summary(ci)), "Wilson")
})

local_adaptive_TS <- function(x, pnull, param) {
  out <- c(
    KS = max(abs(sort(pnull(x)) - (seq_along(x)-0.5)/length(x))),
    MeanAbs = abs(mean(x))
  )
  out
}


test_that("gof_power_adaptive returns structured adaptive output", {
  set.seed(401)

  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(40)
  ralt <- function(mu) rnorm(40, mu)

  out <- gof_power_adaptive(
    pnull, NA, rnull, ralt,
    param_alt=c(0, 0.5),
    TS=local_adaptive_TS,
    B.null=40,
    B.min=40,
    B.max=80,
    B.batch=20,
    target.mcse=0.10,
    maxProcessor=1,
    SuppressMessages=TRUE
  )

  expect_s3_class(out, "Rgof_power_adaptive")
  expect_s3_class(out, "Rgof_power")
  expect_true(all(c(
    "power", "mc.se", "lower", "upper", "rejections",
    "B", "B.null", "critical.values", "target.mcse",
    "converged", "B.min", "B.max", "B.batch"
  ) %in% names(out)))

  expect_true(is.matrix(out$power))
  expect_equal(dim(out$power), dim(out$mc.se))
  expect_equal(dim(out$power), dim(out$lower))
  expect_equal(dim(out$power), dim(out$upper))
  expect_equal(dim(out$power), dim(out$rejections))

  expect_true(out$B >= 40)
  expect_true(out$B <= 80)
  expect_true(all(out$lower <= out$power))
  expect_true(all(out$power <= out$upper))
})


test_that("gof_power_adaptive keeps exact rejection-count identity", {
  set.seed(402)

  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(30)
  ralt <- function(mu) rnorm(30, mu)

  out <- gof_power_adaptive(
    pnull, NA, rnull, ralt,
    param_alt=c(0, 0.8),
    TS=local_adaptive_TS,
    B.null=30,
    B.min=30,
    B.max=60,
    B.batch=15,
    target.mcse=0.12,
    maxProcessor=1,
    SuppressMessages=TRUE
  )

  expect_equal(out$power, out$rejections/out$B)
})


test_that("gof_power_adaptive handles one alternative and retains names", {
  set.seed(403)

  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(30)
  ralt <- function(mu) rnorm(30, mu)

  myTS <- function(x, pnull, param) {
    ans <- abs(mean(x))
    names(ans) <- "MeanAbs"
    ans
  }

  out <- gof_power_adaptive(
    pnull, NA, rnull, ralt,
    param_alt=0.5,
    TS=myTS,
    B.null=30,
    B.min=30,
    B.max=60,
    B.batch=15,
    target.mcse=0.12,
    maxProcessor=1,
    SuppressMessages=TRUE
  )

  expect_false(is.matrix(out$power))
  expect_named(out$power, "MeanAbs")
  expect_named(out$mc.se, "MeanAbs")
  expect_named(out$lower, "MeanAbs")
  expect_named(out$upper, "MeanAbs")
  expect_named(out$rejections, "MeanAbs")
})


test_that("gof_power_adaptive supports a custom p-value routine", {
  set.seed(404)

  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(30)
  ralt <- function(mu) rnorm(30, mu)

  pTS <- function(x, pnull, param) {
    ans <- t.test(x)$p.value
    names(ans) <- "t-test"
    ans
  }

  out <- gof_power_adaptive(
    pnull, NA, rnull, ralt,
    param_alt=0.5,
    TS=pTS,
    With.p.value=TRUE,
    B.min=30,
    B.max=60,
    B.batch=15,
    target.mcse=0.12,
    maxProcessor=1,
    SuppressMessages=TRUE
  )

  expect_identical(out$B.null, 0L)
  expect_null(out$critical.values)
  expect_named(out$power, "t-test")
  expect_equal(out$power, out$rejections/out$B)
})


test_that("gof_power_adaptive validates adaptive simulation controls", {
  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(20)
  ralt <- function(mu) rnorm(20, mu)

  expect_error(
    gof_power_adaptive(
      pnull, NA, rnull, ralt, 0.5,
      B.min=100, B.max=50,
      maxProcessor=1
    ),
    "B.min must not exceed B.max"
  )

  expect_error(
    gof_power_adaptive(
      pnull, NA, rnull, ralt, 0.5,
      B.null=0,
      maxProcessor=1
    ),
    "B.null must be an integer"
  )

  expect_error(
    gof_power_adaptive(
      pnull, NA, rnull, ralt, 0.5,
      target.mcse=0,
      maxProcessor=1
    ),
    "target.mcse"
  )
})


test_that("B.null is independent of the alternative B.max", {
  set.seed(406)

  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(30)
  ralt <- function(mu) rnorm(30, mu)

  out <- gof_power_adaptive(
    pnull, NA, rnull, ralt,
    param_alt=0.5,
    TS=local_adaptive_TS,
    B.null=60,
    B.min=10,
    B.max=20,
    B.batch=5,
    target.mcse=0.20,
    maxProcessor=1,
    SuppressMessages=TRUE
  )

  expect_identical(out$B.null, 60L)
  expect_true(out$B <= 20)
})


test_that("as.data.frame works for adaptive power output", {
  set.seed(405)

  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(25)
  ralt <- function(mu) rnorm(25, mu)

  out <- gof_power_adaptive(
    pnull, NA, rnull, ralt,
    param_alt=c(0, 0.5),
    TS=local_adaptive_TS,
    B.null=20,
    B.min=20,
    B.max=40,
    B.batch=10,
    target.mcse=0.15,
    maxProcessor=1,
    SuppressMessages=TRUE
  )

  d <- as.data.frame(out)

  expect_true(all(c(
    "param_alt", "method", "power", "mc.se",
    "lower", "upper", "rejections"
  ) %in% names(d)))

  expect_equal(nrow(d), length(out$power))
})

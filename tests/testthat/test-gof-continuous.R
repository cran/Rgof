local_continuous_setup <- function(n = 60) {
  list(
    pnull = function(x) pnorm(x),
    rnull = function() rnorm(n),
    ralt = function(m) rnorm(n, m),
    TSextra = list(qnull = function(x) qnorm(x))
  )
}

test_that("gof_test handles continuous data without estimation", {
  set.seed(111)
  z <- local_continuous_setup()
  x <- z$rnull()

  out <- gof_test(x, NA, z$pnull, z$rnull,
                  TSextra = z$TSextra, B = 12, maxProcessor = 1)

  expect_s3_class(out, "Rgof_test")
  expect_true(all(c("results", "statistics", "p.values", "parameters",
                    "metadata", "call") %in% names(out)))
  expect_true(is.data.frame(out$results))
  expect_true(is.numeric(out$statistics))
  expect_true(is.numeric(out$p.values))
  expect_true(length(out$statistics) > 0L)
  expect_identical(names(out$statistics), names(out$p.values))
  expect_true(all(is.finite(out$statistics)))
  expect_true(all(out$p.values >= 0 & out$p.values <= 1))
  expect_equal(out$results$statistic, unname(out$statistics))
  expect_equal(out$results$p.value, unname(out$p.values))
  expect_equal(out$metadata$n, length(x))
  expect_identical(out$metadata$data.type, "continuous")
  expect_true(inherits(out$call, "call"))
})

test_that("gof_test dispatches 3- and 4-argument continuous custom statistics", {
  set.seed(111)
  z <- local_continuous_setup()
  x <- z$rnull()

  myTS1 <- function(x, pnull, param) {
    x <- sort(x)
    Fx <- pnull(x)
    out <- abs(mean(x) - sum(x[-1] * diff(Fx)))
    names(out) <- "M"
    out
  }
  myTS2 <- function(x, pnull, param, TSextra) {
    a <- quantile(x, 1:10 / 11)
    b <- TSextra$qnull(1:10 / 11)
    out <- sum(abs(a - b))
    names(out) <- "A"
    out
  }

  out1 <- gof_test(x, NA, z$pnull, z$rnull, TS = myTS1,
                   B = 12, maxProcessor = 1)
  out2 <- gof_test(x, NA, z$pnull, z$rnull, TS = myTS2,
                   TSextra = z$TSextra, B = 12, maxProcessor = 1)

  expect_s3_class(out1, "Rgof_test")
  expect_s3_class(out2, "Rgof_test")
  expect_named(out1$statistics, "M")
  expect_named(out2$statistics, "A")
  expect_true(out1$p.values >= 0 & out1$p.values <= 1)
  expect_true(out2$p.values >= 0 & out2$p.values <= 1)
})

test_that("gof_test handles continuous parameter estimation", {
  set.seed(111)
  n <- 60
  pnull <- function(x, m) pnorm(x, m)
  rnull <- function(m) rnorm(n, m)
  phat <- function(x) mean(x)
  TSextra <- list(qnull = function(x, m = 0) qnorm(x, m))
  x <- rnorm(n, 1)

  out <- gof_test(x, NA, pnull, rnull, phat = phat,
                  TSextra = TSextra, B = 12, maxProcessor = 1)

  expect_s3_class(out, "Rgof_test")
  expect_identical(names(out$statistics), names(out$p.values))
  expect_true(all(out$p.values >= 0 & out$p.values <= 1))
  expect_length(out$parameters, 1L)
  expect_equal(unname(out$parameters), mean(x))
})

test_that("gof_power handles continuous built-in and custom tests", {
  set.seed(111)
  z <- local_continuous_setup(n = 50)

  base <- gof_power(z$pnull, NA, z$rnull, z$ralt,
                    param_alt = c(0, 0.5), TSextra = z$TSextra,
                    B = 12, maxProcessor = 1)

  myTS <- function(x, pnull, param) {
    out <- abs(mean(x))
    names(out) <- "MeanAbs"
    out
  }
  custom <- gof_power(z$pnull, NA, z$rnull, z$ralt,
                      param_alt = c(0, 0.5), TS = myTS,
                      B = 12, maxProcessor = 1)

  expect_true(is.matrix(base))
  expect_equal(nrow(base), 2L)
  expect_true(all(base >= 0 & base <= 1))
  expect_true(is.matrix(custom))
  expect_equal(dim(custom), c(2L, 1L))
  expect_identical(colnames(custom), "MeanAbs")
  expect_true(all(custom >= 0 & custom <= 1))
})

test_that("gof_power accepts continuous custom p-values", {
  set.seed(111)
  z <- local_continuous_setup(n = 50)
  extra <- list(statistic = FALSE, nbins = 5)

  out <- gof_power(z$pnull, NA, z$rnull, z$ralt,
                   param_alt = c(0, 0.5), TS = myTS_cont,
                   TSextra = extra, With.p.value = TRUE,
                   B = 8, maxProcessor = 1)

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(2L, 1L))
  expect_identical(colnames(out), "Chi_C5")
  expect_true(all(out >= 0 & out <= 1))
})

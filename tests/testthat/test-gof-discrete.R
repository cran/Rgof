local_discrete_setup <- function(n = 200) {
  vals <- 0:20
  list(
    vals = vals,
    pnull = function() pbinom(vals, 20, 0.5),
    rnull = function() table(c(vals, rbinom(n, 20, 0.5))) - 1,
    ralt = function(p) table(c(vals, rbinom(n, 20, p))) - 1,
    TSextra = list(a = 0)
  )
}

test_that("gof_test handles discrete data without estimation", {
  set.seed(111)
  z <- local_discrete_setup()
  x <- z$rnull()

  out <- gof_test(x, z$vals, z$pnull, z$rnull,
                  B = 12, maxProcessor = 1)

  expect_s3_class(out, "Rgof_test")
  expect_true(all(c("results", "statistics", "p.values", "parameters",
                    "metadata", "call") %in% names(out)))
  expect_true(is.data.frame(out$results))
  expect_true(is.numeric(out$statistics))
  expect_true(is.numeric(out$p.values))
  expect_identical(names(out$statistics), names(out$p.values))
  expect_true(all(is.finite(out$statistics)))
  expect_true(all(out$p.values >= 0 & out$p.values <= 1))
  expect_equal(out$results$statistic, unname(out$statistics))
  expect_equal(out$results$p.value, unname(out$p.values))
  expect_equal(out$metadata$n, sum(x))
  expect_identical(out$metadata$data.type, "discrete")
  expect_true(inherits(out$call, "call"))
})

test_that("gof_test dispatches 4- and 5-argument discrete custom statistics", {
  set.seed(111)
  z <- local_discrete_setup()
  x <- z$rnull()

  myTS1 <- function(x, pnull, param, vals) {
    Fx <- pnull()
    out <- abs(sum(x * vals) / sum(x) - sum(vals * diff(c(0, Fx)))) / 20
    names(out) <- "M1"
    out
  }
  myTS2 <- function(x, pnull, param, vals, TSextra) {
    Fx <- pnull()
    out <- abs(sum(x * vals) / sum(x) - sum(vals * diff(c(0, Fx)))) / 20
    names(out) <- "M2"
    out
  }

  out1 <- gof_test(x, z$vals, z$pnull, z$rnull,
                   TS = myTS1, B = 12, maxProcessor = 1)
  out2 <- gof_test(x, z$vals, z$pnull, z$rnull,
                   TS = myTS2, TSextra = z$TSextra,
                   B = 12, maxProcessor = 1)

  expect_s3_class(out1, "Rgof_test")
  expect_s3_class(out2, "Rgof_test")
  expect_true("M1" %in% names(out1$statistics))
  expect_true("M2" %in% names(out2$statistics))
  expect_true(all(out1$p.values >= 0 & out1$p.values <= 1))
  expect_true(all(out2$p.values >= 0 & out2$p.values <= 1))
})

test_that("gof_test handles discrete parameter estimation", {
  set.seed(111)
  n <- 200
  vals <- 0:20
  clamp <- function(p) ifelse(0 < p & p < 1, p, 0.001)
  pnull <- function(p) pbinom(vals, 20, clamp(p))
  rnull <- function(p) table(c(vals, rbinom(n, 20, clamp(p)))) - 1
  phat <- function(x) sum(x * vals) / sum(x) / 20
  x <- table(c(vals, rbinom(n, 20, 0.6))) - 1

  out <- gof_test(x, vals, pnull, rnull, phat = phat,
                  B = 12, maxProcessor = 1)

  expect_s3_class(out, "Rgof_test")
  expect_identical(names(out$statistics), names(out$p.values))
  expect_true(all(out$p.values >= 0 & out$p.values <= 1))
  expect_length(out$parameters, 1L)
  expect_equal(unname(out$parameters), unname(phat(x)))
})

test_that("gof_power handles discrete built-in and custom tests", {
  set.seed(111)
  z <- local_discrete_setup(n = 150)

  base <- gof_power(z$pnull, z$vals, z$rnull, z$ralt,
                    param_alt = c(0.5, 0.6),
                    B = 12, maxProcessor = 1)

  myTS <- function(x, pnull, param, vals, TSextra) {
    Fx <- pnull()
    out <- abs(sum(x * vals) / sum(x) - sum(vals * diff(c(0, Fx)))) / 20
    names(out) <- "M2"
    out
  }
  custom <- gof_power(z$pnull, z$vals, z$rnull, z$ralt,
                      param_alt = c(0.5, 0.6), TS = myTS,
                      TSextra = z$TSextra,
                      B = 12, maxProcessor = 1)

  expect_true(is.matrix(base))
  expect_equal(nrow(base), 2L)
  expect_true(all(base >= 0 & base <= 1))
  expect_true(is.matrix(custom))
  expect_equal(dim(custom), c(2L, 1L))
  expect_identical(colnames(custom), "M2")
  expect_true(all(custom >= 0 & custom <= 1))
})

test_that("gof_power accepts discrete custom p-values", {
  set.seed(111)
  z <- local_discrete_setup(n = 150)
  extra <- list(statistic = FALSE, nbins = 5)

  out <- gof_power(z$pnull, z$vals, z$rnull, z$ralt,
                   param_alt = c(0.5, 0.6), TS = myTS_disc,
                   TSextra = extra, With.p.value = TRUE,
                   B = 8, maxProcessor = 1)

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(2L, 1L))
  expect_identical(colnames(out), "Chi_D")
  expect_true(all(out >= 0 & out <= 1))
})

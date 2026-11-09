test_that("hybrid_test runs for continuous and discrete data", {
  skip_if_not_installed("R2sample", minimum_version="5.0.1")
  set.seed(111)
  rnull <- function() rnorm(100)
  out <- hybrid_test(rnull(), rnull=rnull, B=20, maxProcessor=1,
                     SuppressMessages=TRUE)
  expect_true(all(out$p.values >= 0 & out$p.values <= 1))
  # composite null hypothesis, data clearly not normal
  rnull <- function(p) rnorm(100, p[1], p[2])
  phat <- function(x) c(mean(x), sd(x))
  out <- hybrid_test(rt(100, 1), rnull=rnull, phat=phat, B=20, maxProcessor=1,
                     SuppressMessages=TRUE, doMethods=c("KS", "AD"))
  expect_true(all(out$p.values < 0.1))
  # discrete data with parameter estimation
  vals <- 0:10
  rnull <- function(p) table(c(0:10, rbinom(1000, 10, p))) - 1
  phat <- function(x) mean(rep(0:10, x))/10
  out <- hybrid_test(rnull(0.5), vals, rnull, phat, B=20, maxProcessor=1,
                     SuppressMessages=TRUE)
  expect_true(all(out$p.values >= 0 & out$p.values <= 1))
})

test_that("hybrid bootstrap re-estimates the parameters", {
  # the Monte Carlo data of each bootstrap sample has to be generated with
  # the parameters estimated from the bootstrap data, not from the observed data
  rnull <- function(p) rnorm(50, p[1], 1)
  phat <- function(x) mean(x)
  hyb <- Rgof:::make_hybrid_generators(rnorm(50, 10), NA, rnull, phat, 1, TRUE)
  set.seed(1)
  b <- hyb$rnullT(list(x=rnorm(50, 10)))
  expect_lt(abs(mean(b$y) - mean(b$x)), 1)
  expect_equal(length(b$y), 50L)
})

test_that("hybrid_power runs", {
  skip_if_not_installed("R2sample", minimum_version="5.0.1")
  set.seed(112)
  pwr <- hybrid_power(function() rnorm(100), function(mu) rnorm(100, mu),
                      c(0, 1), B=20, maxProcessor=1, SuppressMessages=TRUE,
                      doMethods=c("KS", "AD"))
  expect_true(all(pwr >= 0 & pwr <= 1))
  expect_equal(ncol(pwr), 2L)
})

test_that("summary of hybrid_power gives a study summary", {
  skip_if_not_installed("R2sample", minimum_version="5.0.1")
  set.seed(113)
  pwr <- hybrid_power(function() rnorm(100), function(mu) rnorm(100, mu),
                      c(0, 1), nMC=2, B=20, maxProcessor=1, SuppressMessages=TRUE,
                      doMethods=c("KS", "AD"))
  expect_true(is.matrix(pwr))
  s <- summary(pwr)
  expect_s3_class(s, "summary.Rgof_power")
  expect_equal(s$metadata$nMC, 200)
  expect_output(print(s), "Rgof hybrid power study summary")
  expect_output(print(s), "Monte Carlo sample size: 200")
})

test_that("hybrid_test returns an Rgof_test object with a summary", {
  skip_if_not_installed("R2sample", minimum_version="5.0.1")
  set.seed(114)
  vals <- 0:10
  rnull <- function(p) table(c(0:10, rbinom(1000, 10, p))) - 1
  phat <- function(x) mean(rep(0:10, x))/10
  out <- hybrid_test(rnull(0.5), vals, rnull, phat, B=20, maxProcessor=1,
                     SuppressMessages=TRUE)
  expect_s3_class(out, "Rgof_test")
  expect_true(is.data.frame(out$results))
  expect_equal(out$results$p.value, unname(out$p.values))
  expect_true(out$metadata$hybrid)
  expect_equal(out$metadata$nMC, 1000)
  expect_output(print(out), "Rgof hybrid goodness-of-fit tests")
  expect_output(print(summary(out)), "Monte Carlo sample size: 1000")
  expect_output(print(summary(out)), "Estimated parameters")
})

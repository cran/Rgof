structured_continuous_result <- function(B = 10) {
  n <- 60
  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(n)
  TSextra <- list(qnull = function(x) qnorm(x))

  x <- rnull()
  gof_test(x, NA, pnull, rnull,
           TSextra = TSextra, B = B, maxProcessor = 1)
}

test_that("Rgof_test has the expected structured components", {
  set.seed(201)
  out <- structured_continuous_result()

  expect_s3_class(out, "Rgof_test")
  expect_true(all(c("results", "statistics", "p.values", "parameters",
                    "metadata", "call") %in% names(out)))

  expect_s3_class(out$results, "data.frame")
  expect_identical(names(out$results), c("method", "statistic", "p.value"))
  expect_equal(nrow(out$results), length(out$statistics))
  expect_equal(nrow(out$results), length(out$p.values))

  expect_equal(out$results$method, names(out$statistics))
  expect_equal(out$results$statistic, unname(out$statistics))
  expect_equal(out$results$p.value, unname(out$p.values))
})

test_that("Rgof_test preserves legacy statistics and p.values access", {
  set.seed(202)
  out <- structured_continuous_result()

  expect_true(is.numeric(out$statistics))
  expect_true(is.numeric(out$p.values))
  expect_identical(names(out$statistics), names(out$p.values))

  expect_equal(
    out$statistics,
    setNames(out$results$statistic, out$results$method)
  )
  expect_equal(
    out$p.values,
    setNames(out$results$p.value, out$results$method)
  )
})

test_that("Rgof_test metadata records the calculation", {
  set.seed(203)
  out <- structured_continuous_result(B = 11)

  expect_identical(out$metadata$data.type, "continuous")
  expect_equal(out$metadata$n, 60)
  expect_equal(out$metadata$B.requested, 11)
  expect_true(out$metadata$B.used >= out$metadata$B.requested)
  expect_equal(out$metadata$B.used %% out$metadata$maxProcessor, 0)
  expect_true(is.numeric(out$metadata$typeTS))
  expect_true(is.numeric(out$metadata$nbins))
  expect_true(inherits(out$call, "call"))
})

test_that("as.data.frame returns the results table", {
  set.seed(204)
  out <- structured_continuous_result()
  tab <- as.data.frame(out)

  expect_s3_class(tab, "data.frame")
  expect_identical(tab, out$results)
})

test_that("summary returns an Rgof_test summary object", {
  set.seed(205)
  out <- structured_continuous_result()
  s <- summary(out)

  expect_s3_class(s, "summary.Rgof_test")
  expect_true(all(c("call", "results", "parameters", "metadata") %in% names(s)))
  expect_identical(s$results, out$results)
  expect_identical(s$parameters, out$parameters)
  expect_identical(s$metadata, out$metadata)
  expect_identical(s$call, out$call)
})

test_that("print methods return their input invisibly", {
  set.seed(206)
  out <- structured_continuous_result()
  s <- summary(out)

  printed_out <- capture.output(ret1 <- print(out))
  printed_summary <- capture.output(ret2 <- print(s))

  expect_gt(length(printed_out), 0L)
  expect_gt(length(printed_summary), 0L)
  expect_identical(ret1, out)
  expect_identical(ret2, s)
})

test_that("gof_test returns structured output when B is zero", {
  set.seed(207)
  out <- structured_continuous_result(B = 0)

  expect_s3_class(out, "Rgof_test")
  expect_s3_class(out$results, "data.frame")
  expect_true(all(is.na(out$p.values)))
  expect_true(all(is.na(out$results$p.value)))
  expect_equal(out$metadata$B.requested, 0)
  expect_equal(out$metadata$B.used, 0)
  expect_true(is.na(out$metadata$maxProcessor))
})

test_that("no-estimation and estimation cases record parameters correctly", {
  set.seed(208)

  no_est <- structured_continuous_result()
  expect_length(no_est$parameters, 0L)

  n <- 60
  pnull <- function(x, m) pnorm(x, m)
  rnull <- function(m) rnorm(n, m)
  phat <- function(x) mean(x)
  x <- rnorm(n, 1)

  est <- gof_test(x, NA, pnull, rnull, phat = phat,
                  B = 10, maxProcessor = 1)

  expect_s3_class(est, "Rgof_test")
  expect_length(est$parameters, 1L)
  expect_equal(unname(est$parameters), mean(x))
})

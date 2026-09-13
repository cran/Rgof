test_that("case_studies lists 60 unique case studies", {
  studies <- case_studies(ReturnCaseNames = TRUE)

  expect_type(studies, "character")
  expect_length(studies, 60L)
  expect_length(unique(studies), 60L)
  expect_identical(studies[c(1, 21, 41, 60)],
                   c("uniform", "uniform-linear", "normal.t.est", 
                     "gammamixture2.est"))
})

test_that("case_studies preserves parameter-estimation signatures", {
  no_est <- case_studies("normal.t", "mixed", "cont", n = 40)
  est <- case_studies("normal.t.est", "mixed", "cont", n = 40)

  expect_length(formals(no_est$pnull), 1L)
  expect_length(formals(no_est$dnull), 1L)
  expect_length(formals(no_est$rnull), 0L)

  expect_length(formals(est$pnull), 2L)
  expect_length(formals(est$dnull), 2L)
  expect_length(formals(est$rnull), 1L)
  expect_true(is.function(est$phat))
})

test_that("case_studies returns usable continuous and discrete studies", {
  set.seed(111)

  cont <- case_studies("uniform", "mean", "cont", n = 40)
  expect_true(all(c("pnull", "rnull", "ralt", "phat", "dnull", "Range") %in%
                    names(cont)))
  expect_length(cont$rnull(), 40L)
  expect_length(cont$ralt(), 40L)
  expect_true(all(cont$pnull(c(0, 0.5, 1)) >= 0 &
                  cont$pnull(c(0, 0.5, 1)) <= 1))

  disc <- case_studies("uniform", "mean", "disc", n = 40, nbins = 10)
  expect_true(all(c("pnull", "rnull", "ralt", "phat", "vals") %in%
                    names(disc)))
  expect_length(disc$vals, 10L)
  expect_length(disc$rnull(), 10L)
  expect_length(disc$ralt(), 10L)
})

test_that("pareto case exposes qnull by name", {
  pareto <- case_studies("pareto.est", "mixed", "cont", n = 40)

  expect_true("qnull" %in% names(pareto$TSextra))
  expect_true(is.function(pareto$TSextra$qnull))
})

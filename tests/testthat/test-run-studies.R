test_that("run.studies returns stored continuous and discrete comparisons", {
  set.seed(111)

  cont <- NULL
  expect_message(
    invisible(capture.output(
      cont <- run.studies(studies = 1, types = 1, data_type = "cont", B = 5)
    )),
    "Average rank of a method"
  )

  disc <- NULL
  expect_message(
    invisible(capture.output(
      disc <- run.studies(studies = 1, types = 1, data_type = "disc", B = 5)
    )),
    "Average rank of a method"
  )

  expect_true(is.matrix(cont))
  expect_true(is.matrix(disc))
  expect_equal(nrow(cont), 1L)
  expect_equal(nrow(disc), 1L)
  expect_true(ncol(cont) > 1L)
  expect_true(ncol(disc) > 1L)
})

test_that("run.studies accepts custom continuous and discrete statistics", {
  set.seed(111)

  cont_builtin <- NULL
  expect_message(
    invisible(capture.output(
      cont_builtin <- run.studies(
        studies = 21, types = 5, data_type = "cont",
        n = 80, B = 5
      )
    )),
    "Average rank of a method"
  )

  disc_builtin <- NULL
  expect_message(
    invisible(capture.output(
      disc_builtin <- run.studies(
        studies = 21, types = 5, data_type = "disc",
        n = 80, B = 5
      )
    )),
    "Average rank of a method"
  )

  cont_extra <- list(statistic = TRUE, nbins = 5)
  cont <- NULL
  expect_message(
    invisible(capture.output(
      cont <- run.studies(
        studies = 21, types = 5, data_type = "cont",
        TS = myTS_cont, TSextra = cont_extra,
        n = 80, B = 5
      )
    )),
    "Average rank of a method"
  )

  disc_extra <- list(statistic = TRUE, nbins = 5)
  disc <- NULL
  expect_message(
    invisible(capture.output(
      disc <- run.studies(
        studies = 21, types = 5, data_type = "disc",
        TS = myTS_disc, TSextra = disc_extra,
        n = 80, B = 5
      )
    )),
    "Average rank of a method"
  )

  expect_true(is.matrix(cont_builtin))
  expect_true(is.matrix(disc_builtin))
  expect_true(is.matrix(cont))
  expect_true(is.matrix(disc))

  expect_equal(nrow(cont), nrow(cont_builtin))
  expect_equal(nrow(disc), nrow(disc_builtin))

  # A user-supplied statistic is prepended to the built-in statistics.
  # run.studies currently does not preserve the custom statistic's column
  # name because it obtains names(pwr) from a matrix returned by gof_power().
  expect_gt(ncol(cont), ncol(cont_builtin))
  expect_gt(ncol(disc), ncol(disc_builtin))
})

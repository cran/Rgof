# Checks of the built-in test statistics against their definitions.

test_that("continuous statistics agree with published formulas", {
  set.seed(1)
  for (r in 1:30) {
    n <- sample(5:80, 1)
    x <- sort(rnorm(n, 0.2))
    u <- pnorm(x); i <- 1:n
    Dp <- max(i/n - u); Dm <- max(u - (i - 1)/n)
    W2 <- 1/(12*n) + sum((u - (2*i - 1)/(2*n))^2)
    a <- Rgof:::TS_cont(x, function(x) pnorm(x), -99, function(x) qnorm(x))
    expect_equal(unname(a["KS"]), max(Dp, Dm), tolerance = 1e-10)
    expect_equal(unname(a["K"]), Dp + Dm, tolerance = 1e-10)
    expect_equal(unname(a["CvM"]), W2, tolerance = 1e-10)
    expect_equal(unname(a["W"]), W2 - n*(mean(u) - 0.5)^2, tolerance = 1e-10)
    expect_equal(unname(a["AD"]),
                 -n - mean((2*i - 1)*(log(u) + log(1 - rev(u)))), tolerance = 1e-10)
    w1 <- sum(sapply(1:n, function(k)
      integrate(function(t) abs(x[k] - qnorm(t)), (k - 1)/n, k/n,
                subdivisions = 1000)$value))
    expect_equal(unname(a["Wassp1"]), w1, tolerance = 1e-3)
  }
})

test_that("discrete statistics agree with their definitions", {
  set.seed(2)
  vals <- 0:10
  Fv <- pbinom(vals, 10, 0.5); p <- diff(c(0, Fv))
  t <- (p + c(p[-1], p[1]))/2   # Choulakian, Lockhart and Stephens (1994)
  k <- length(vals)
  for (r in 1:20) {
    n <- sample(20:300, 1)
    x <- as.integer(table(factor(rbinom(n, 10, 0.53), levels = vals)))
    Z <- cumsum(x)/n - Fv
    a <- Rgof:::TS_disc(x, function() Fv, -99, vals)
    expect_equal(unname(a["KS"]), max(abs(Z)), tolerance = 1e-10)
    expect_equal(unname(a["K"]), max(Z, 0) + max(-Z, 0), tolerance = 1e-10)
    expect_equal(unname(a["CvM"]), n*sum(Z^2*t), tolerance = 1e-10)
    expect_equal(unname(a["W"]), n*sum((Z - sum(Z*t))^2*t), tolerance = 1e-10)
    expect_equal(unname(a["AD"]), n*sum((Z^2*t/(Fv*(1 - Fv)))[-k]), tolerance = 1e-10)
  }
  # agreement with the implementation in package dgof
  skip_if_not_installed("dgof")
  x <- as.integer(table(factor(rbinom(200, 10, 0.53), levels = vals)))
  a <- Rgof:::TS_disc(x, function() Fv, -99, vals)
  xx <- rep(vals, x)
  sf <- stepfun(vals, c(0, Fv))
  expect_equal(unname(a["KS"]), unname(dgof::ks.test(xx, sf)$statistic), tolerance = 1e-10)
  expect_equal(unname(a["CvM"]), unname(dgof::cvm.test(xx, sf, type = "W2")$statistic), tolerance = 1e-8)
  expect_equal(unname(a["W"]), unname(dgof::cvm.test(xx, sf, type = "U2")$statistic), tolerance = 1e-8)
  expect_equal(unname(a["AD"]), unname(dgof::cvm.test(xx, sf, type = "A2")$statistic), tolerance = 1e-8)
})

test_that("discrete smooth tests and Berk-Jones agree with their definitions", {
  # reference implementation: polynomials orthonormal w.r.t. p via qr()
  ref <- function(x, Fv) {
    k <- length(x); n <- sum(x)
    p <- pmax(diff(c(0, Fv)), 0); u <- c(0, Fv[-k]) + p/2
    J <- p > 0; K <- min(10, sum(J) - 1); t <- 2*u[J] - 1
    P <- matrix(1, sum(J), K + 1); P[, 2] <- t
    if (K >= 2) for (j in 2:K) P[, j + 1] <- ((2*j - 1)*t*P[, j] - (j - 1)*P[, j - 1])/j
    g <- qr.Q(qr(sqrt(p[J])*P))/sqrt(p[J])
    S <- cumsum((colSums(x[J]*g)[-1]/sqrt(n))^2)
    I <- which(Fv[-k] > 0 & Fv[-k] < 1)
    a <- (cumsum(x)/n)[I]; f <- Fv[I]
    kl <- ifelse(a == 0, -log(1 - f), ifelse(a == 1, -log(f),
           a*log(a/f) + (1 - a)*log((1 - a)/(1 - f))))
    c("NS-B" = S[which.max(S - (1:K)*log(n))], "NS-A" = S[which.max(S - 2*(1:K))],
      BJ = n*max(kl))
  }
  set.seed(5)
  vals <- 0:10
  Fv <- pbinom(vals, 10, 0.5)
  for (r in 1:20) {
    n <- sample(20:300, 1)
    x <- as.integer(table(factor(rbinom(n, 10, 0.53), levels = vals)))
    a <- Rgof:::TS_disc(x, function() Fv, -99, vals)
    expect_equal(a[c("NS-B", "NS-A", "BJ")], ref(x, Fv), tolerance = 1e-8)
  }
  # a value with null probability 0, and only three values
  x <- c(0L, 5L, 12L, 30L, 20L, 8L, 5L); Fv <- c(0, 0.1, 0.3, 0.6, 0.85, 0.95, 1)
  expect_equal(Rgof:::TS_disc(x, function() Fv, -99, 1:7)[c("NS-B", "NS-A", "BJ")],
               ref(x, Fv), tolerance = 1e-8)
  expect_equal(unname(Rgof:::TS_disc(c(25L, 55L, 20L), function() c(0.3, 0.8, 1), -99,
                      1:3)[c("NS-B", "NS-A", "BJ")]),
               c(0.5714285714, 0.5714285714, 0.6164264417), tolerance = 1e-8)
})

test_that("weighted statistics with equal weights equal unweighted ones", {
  set.seed(3)
  x <- rnorm(40)
  a <- Rgof:::TSw_cont(x, function(x) pnorm(x), -99, rep(1, 40))
  b <- Rgof:::TS_cont(x, function(x) pnorm(x), -99, function(x) qnorm(x))
  expect_equal(unname(a[c("KS", "K", "CvM", "AD")]),
               unname(b[c("KS", "K", "CvM", "AD")]), tolerance = 1e-8)
})

test_that("adjusted p value uses the observed data", {
  set.seed(4)
  pnull <- function(x) pnorm(x)
  rnull <- function() rnorm(100)
  x <- rnorm(100, 1)
  out <- suppressMessages(
    gof_test_adjusted_pvalue(x, NA, pnull, rnull, B = c(200, 100),
                             maxProcessor = 1, doMethods = c("KS", "AD"))
  )
  expect_type(out, "list")
  expect_true(all(out$p.values == 0))
  # with the observed data all p values are 0; the adjusted p value is small but
  # can be positive, because simulated data sets can also have a smallest p value of 0
  expect_lt(out$adjusted.p.value, 0.1)
})

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{power_disc}
\alias{power_disc}
\title{find power of gof tests for discrete data}
\usage{
power_disc(
  pnull,
  rnull,
  vals,
  ralt,
  param_alt,
  phat,
  TS,
  nbins = as.integer(c(100, 10)),
  rate = 0,
  B = as.integer(c(1000, 1000)),
  alpha = 0.05,
  minexpcount = 2
)
}
\arguments{
\item{pnull}{R function (cdf)}

\item{rnull}{R function (generate data under null hypothesis)}

\item{vals}{vector of values of discrete random variable}

\item{ralt}{R function to generate data under alternative}

\item{param_alt}{parameters of function ralt}

\item{phat}{function to estimate parameters from the data}

\item{TS}{function to calculate test statistics}

\item{nbins}{=c(100,10) number of bins to use}

\item{rate}{rate of Poisson if sample size is random, 0 otherwise}

\item{B}{=c(1000, 1000) Number of simulation runs for power and null distribution}

\item{alpha}{=0.05, type I error of test}

\item{minexpcount}{=2 minimal expected bin count required}
}
\value{
A matrix of powers
}
\description{
find power of gof tests for discrete data
}
\keyword{internal}

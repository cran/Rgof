## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error = TRUE
)

## ----setup--------------------------------------------------------------------
library(Rgof)

## -----------------------------------------------------------------------------
set.seed(123)

## -----------------------------------------------------------------------------
vals=0:20 #possible values of random variable
pnull=function()  pbinom(0:20, 20, 0.5)  # cumulative distribution function (cdf)
rnull = function() table(c(0:20, rbinom(1000, 20, 0.5)))-1 
# generate data under the null hypothesis, make sure that vector of counts has same length as vals, possibly 0.

## -----------------------------------------------------------------------------
x = rnull()
gof_test_disc(x, pnull, rnull, vals, B=1000, doMethod = "all")

## ----d2-----------------------------------------------------------------------
x = table(c(0:20, rbinom(1000, 20, 0.55)))-1
#true p is 0.55, not 0.5
gof_test_disc(x, pnull, rnull, vals, B=1000, doMethod = "all")$p.value

## ----r1, eval=FALSE-----------------------------------------------------------
#  rnull = function() table(c(0:20, rbinom(rpois(1, 650), 20, 0.5)))-1
#  x = rnull()
#  gof_test_disc(x, pnull, rnull, vals, rate=650, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
vals=0:20
pnull=function(p=0.5)  pbinom(0:20, 20, ifelse(p>0&&p<1, p, 0.5))  
rnull = function(p=0.5) table(c(0:20, rbinom(1000, 20, p)))-1
phat = function(x) sum(0:20*x)/sum(x)/20

## ----ed1----------------------------------------------------------------------
x = table(c(0:20, rbinom(1000, 20, 0.5)))-1  
gof_test_disc(x, pnull, rnull, vals, phat=phat, B=1000, doMethod = "all")$p.value

## ----d3-----------------------------------------------------------------------
x = table(c(0:20, rbinom(1000, 20, 0.55)))-1 
# p is not 0.5, but data is still from a binomial distribution with n=20
gof_test_disc(x, pnull, rnull, vals, phat=phat, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
x = table(c(rep(0:20, 5), rbinom(1000-21*5, 20, 0.53))) 
# data has to many small and large values to be from a binomial
gof_test_disc(x, pnull, rnull, vals, phat=phat, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
rnull = function() {
  y = rexp(2500, 1) # Exp(1) data
  y = y[y<2][1:1500] # 1500 events on 0-2
  bins = 0:40/20 # binning
  hist(y, bins, plot=FALSE)$counts # find bin counts
}
x = rnull()
bins = 0:40/20
vals = (bins[-1]+bins[-21])/2
pnull = function() {
   bins = 1:40/20
   pexp(bins, 1)/pexp(2, 1)
}
  
gof_test_disc(x, pnull, rnull, vals, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
pnull = function(x) pnorm(x)
qnull = function(x) qnorm(x)
rnull = function()  rnorm(1000)

## -----------------------------------------------------------------------------
x = rnorm(1000)
gof_test_cont(x, pnull, rnull, qnull, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5) 
gof_test_cont(x, pnull, rnull, qnull, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
pnull = function(x, p=0) pnorm(x, p)
qnull = function(x, p=0) qnorm(x, p)
rnull = function(p)  rnorm(1000, p)
phat = function(x) mean(x)

## -----------------------------------------------------------------------------
x = rnorm(1000) 
gof_test_cont(x, pnull, rnull, qnull, phat=phat, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5) 
gof_test_cont(x, pnull, rnull, qnull, phat=phat)$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5, 2) 
gof_test_cont(x, pnull, rnull, qnull, phat=phat, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
pnull = function(x, p=c(0, 1)) pnorm(x, p[1], ifelse(p[2]>0, p[2], 0.001))
qnull = function(x, p=c(0, 1)) qnorm(x, p[1], ifelse(p[2]>0, p[2], 0.001))
rnull = function(p=c(0, 1))  rnorm(1000, p[1], ifelse(p[2]>0, p[2], 0.001))
phat = function(x) c(mean(x), sd(x))

## -----------------------------------------------------------------------------
x = rnorm(1000) 
gof_test_cont(x, pnull, rnull, qnull, phat=phat, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5) 
gof_test_cont(x, pnull, rnull, qnull, phat=phat, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5, 2) 
gof_test_cont(x, pnull, rnull, qnull, phat=phat, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
x = rt(1000, 2) 
gof_test_cont(x, pnull, rnull, qnull, phat=phat, B=1000, doMethod = "all")$p.value

## -----------------------------------------------------------------------------
vals = 0:10
pnull = function() pbinom(0:10, 10, 0.5)
rnull =function () table(c(0:10, rbinom(100, 10, 0.5)))-1
ralt =function (p) table(c(0:10, rbinom(100, 10, p)))-1
P=gof_power_disc(pnull, rnull, vals, ralt, 
  param_alt=seq(0.5, 0.6, 0.01),  B=c(500, 500), nbins=c(11, 5))
plot_power(P, "p", Smooth=FALSE)

## -----------------------------------------------------------------------------
vals = 0:10
pnull = function(p) pbinom(0:10, 10, p)
rnull = function (p) table(c(0:10, rbinom(100, 10, p)))-1
phat = function(x) sum(0:10*x)/1000

## -----------------------------------------------------------------------------
ralt =function (p) table(c(0:10, rbinom(100, 10, p)))-1
gof_power_disc(pnull, rnull, vals, ralt, c(0.5, 0.6), phat,
        B=c(100, 200), nbins=c(11, 5), maxProcessors = 2)


## -----------------------------------------------------------------------------
ralt =function (p) table(c(rep(0:10, 2), rbinom(100, 10, p)))
gof_power_disc(pnull, rnull, vals, ralt, 0.5, phat,
        B=c(100, 200), nbins=c(11, 5), maxProcessors = 2)

## -----------------------------------------------------------------------------
pnull = function(x) pnorm(x)
qnull = function(x) qnorm(x)
rnull = function() rnorm(100)
ralt = function(mu=0) rnorm(100, mu)
gof_power_cont(pnull, rnull, qnull, ralt, c(0, 1), B=c(100, 200))

## -----------------------------------------------------------------------------
pnull = function(x, p=c(0,1)) pnorm(x, p[1], ifelse(p[2]>0, p[2], 0.01))
qnull = function(x, p=c(0,1)) qnorm(x, p[1], ifelse(p[2]>0, p[2], 0.01))
rnull = function(p=c(0,1)) rnorm(500, p[1], p[2])
ralt = function(mu=0) rnorm(100, mu)
phat = function(x) c(mean(x), sd(x))
gof_power_cont(pnull, rnull, qnull, ralt, c(0, 1), phat, B=c(100, 200),
               maxProcessor=2)

## -----------------------------------------------------------------------------
ralt = function(df=1) {
# t distribution truncated at +- 5  
  x=rt(1000, df)
  x=x[abs(x)<5]
  x[1:100]
}  
gof_power_cont(pnull, rnull, qnull, ralt, c(2, 50), phat, Range=c(-5,5), B=c(100, 200), maxProcessor=2)

## -----------------------------------------------------------------------------
myTS = function(x, Fx, param, qnull) {
   x=sort(x)
   Fx=Fx= (1:length(x))/length(x)
   out = sum(abs(x-Fx))
   names(out) = "CvM alt"
   out
}

## -----------------------------------------------------------------------------
pnull = function(x) punif(x)
qnull = function(x) qunif(x)
rnull = function() runif(500)
x = rnull()
Rgof::gof_test_cont(x, pnull, rnull, qnull, TS=myTS)

## -----------------------------------------------------------------------------
ralt = function(slope) {
  if(slope==0) y=runif(500)
    else y=(slope-1+sqrt((1-slope)^2+4*slope* runif(500)))/2/slope
}

## -----------------------------------------------------------------------------
gof_power_cont(pnull, rnull, qnull, ralt, TS=myTS, param_alt=seq(0, 0.5, length=5), Range=c(0,1))

## -----------------------------------------------------------------------------
myTS = function(x, p, nm=0, vals) {
   z= sum(x*vals)/sum(x) - sum(p*vals)
   names(z) = "Means"
   z
}

## -----------------------------------------------------------------------------
vals=0:10
pnull = function() pbinom(0:10, 10, 0.5)
rnull = function() table(c(0:10, rbinom(1000, 10, 0.5)))-1
x = rnull()
gof_test_disc(x, pnull, rnull, vals, TS=myTS)

## -----------------------------------------------------------------------------
ralt = function(p) table(c(0:10, rbinom(1000, 10, p)))-1
gof_power_disc(pnull, rnull, vals, ralt, TS=myTS, param_alt=c(0.5, 0.51, 0.52, 0.53))


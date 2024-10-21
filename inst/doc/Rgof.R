## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error = TRUE
)

## ----setup--------------------------------------------------------------------
library(Rgof)
Bsim = c(100, 200) #Number of Simulation Runs

## -----------------------------------------------------------------------------
set.seed(123)

## -----------------------------------------------------------------------------
vals=0:20 #possible values of random variable
pnull=function()  pbinom(0:20, 20, 0.5)  # cumulative distribution function (cdf)
rnull = function() table(c(0:20, rbinom(1000, 20, 0.5)))-1 
# generate data under the null hypothesis, make sure that vector of counts has same length as vals, possibly 0.

## -----------------------------------------------------------------------------
x = rnull()
# Basic Test
gof_test(x, vals, pnull, rnull, B=1000)
#Test with adjusted overall p value
gof_test_adjusted_pvalue(x, vals, pnull, rnull, B=c(1000, 500))

## ----d2-----------------------------------------------------------------------
x = table(c(0:20, rbinom(1000, 20, 0.55)))-1
#true p is 0.55, not 0.5
# Basic Test
gof_test(x, vals, pnull, rnull, B=1000, doMethod = "all")$p.value
#Test with adjusted overall p value
gof_test_adjusted_pvalue(x, vals, pnull, rnull, B=c(1000, 500))

## ----r1, eval=FALSE-----------------------------------------------------------
#  rnull = function() table(c(0:20, rbinom(rpois(1, 650), 20, 0.5)))-1
#  x = rnull()
#  gof_test(x, vals, pnull, rnull, rate=650, B=1000)$p.value

## -----------------------------------------------------------------------------
vals=0:20
pnull=function(p=0.5)  pbinom(0:20, 20, ifelse(p>0&&p<1, p, 0.5))  
rnull = function(p=0.5) table(c(0:20, rbinom(1000, 20, p)))-1
phat = function(x) sum(0:20*x)/sum(x)/20

## ----ed1----------------------------------------------------------------------
x = table(c(0:20, rbinom(1000, 20, 0.5)))-1  
gof_test(x, vals, pnull, rnull, phat=phat, B=1000)$p.value

## ----d3-----------------------------------------------------------------------
x = table(c(0:20, rbinom(1000, 20, 0.55)))-1 
# p is not 0.5, but data is still from a binomial distribution with n=20
gof_test(x, vals, pnull, rnull, phat=phat, B=1000)$p.value

## -----------------------------------------------------------------------------
x = table(c(rep(0:20, 5), rbinom(1000-21*5, 20, 0.53))) 
# data has to many small and large values to be from a binomial
gof_test(x, vals, pnull, rnull, phat=phat, B=1000)$p.value

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
  
gof_test(x, vals, pnull, rnull)$p.value

## -----------------------------------------------------------------------------
pnull = function(x) pnorm(x)
rnull = function()  rnorm(1000)
TSextra = list(qnull=function(x) qnorm(x)) #optional quantile function used by chi square tests and Wassp1 test.

## -----------------------------------------------------------------------------
x = rnorm(1000)
#Basic Tests
gof_test(x, NA, pnull, rnull, B=1000, TSextra=TSextra)$p.value
#Adjusted p value
gof_test_adjusted_pvalue(x, NA, pnull, rnull, B=c(1000,500), TSextra=TSextra)

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5) 
gof_test(x, NA, pnull, rnull, B=1000, TSextra=TSextra)$p.value

## -----------------------------------------------------------------------------
pnull = function(x, p=0) pnorm(x, p)
TSextra = list(qnull = function(x, p=0) qnorm(x, p))
rnull = function(p)  rnorm(1000, p)
phat = function(x) mean(x)

## -----------------------------------------------------------------------------
x = rnorm(1000) 
gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra, B=1000)$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5) 
gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra)$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5, 2) 
gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra, B=1000)$p.value

## -----------------------------------------------------------------------------
pnull = function(x, p=c(0, 1)) pnorm(x, p[1], ifelse(p[2]>0, p[2], 0.001))
TSextra = list(qnull = function(x, p=c(0, 1)) qnorm(x, p[1], ifelse(p[2]>0, p[2], 0.001)))
rnull = function(p=c(0, 1))  rnorm(1000, p[1], ifelse(p[2]>0, p[2], 0.001))
phat = function(x) c(mean(x), sd(x))

## -----------------------------------------------------------------------------
x = rnorm(1000) 
gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra, B=1000)$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5) 
gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra, B=1000)$p.value

## -----------------------------------------------------------------------------
x = rnorm(1000, 0.5, 2) 
gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra, B=1000)$p.value

## -----------------------------------------------------------------------------
x = rt(1000, 2) 
gof_test(x, NA, pnull, rnull, phat=phat, TSextra=TSextra, B=1000)$p.value

## ----pow1---------------------------------------------------------------------
vals = 0:10
pnull = function() pbinom(0:10, 10, 0.5)
rnull =function () table(c(0:10, rbinom(100, 10, 0.5)))-1
ralt =function (p=0.5) table(c(0:10, rbinom(100, 10, p)))-1
P=gof_power(pnull, vals, rnull, ralt, 
  param_alt=seq(0.5, 0.6, 0.02),  B=Bsim, nbins=c(11, 5))
plot_power(P, "p", Smooth=FALSE)

## -----------------------------------------------------------------------------
vals = 0:10
pnull = function(p=0.5) pbinom(0:10, 10, ifelse(0<p&p<1,p,0.001))
rnull = function (p=0.5) table(c(0:10, rbinom(100, 10, ifelse(0<p&p<1,p,0.001))))-1
phat = function(x) sum(0:10*x)/1000

## ----pow2---------------------------------------------------------------------
ralt =function (p=0.5) table(c(0:10, rbinom(100, 10, p)))-1
gof_power(pnull, vals, rnull, ralt, c(0.5, 0.6), phat=phat,
        B=Bsim, nbins=c(11, 5), maxProcessors = 2)


## -----------------------------------------------------------------------------
ralt =function (p=0.5) table(c(rep(0:10, 2), rbinom(100, 10, p)))
gof_power(pnull, vals, rnull, ralt, 0.5, phat=phat,
        B=Bsim, nbins=c(11, 5), maxProcessors = 2)

## -----------------------------------------------------------------------------
pnull = function(x) pnorm(x)
TSextra = list(qnull = function(x) qnorm(x))
rnull = function() rnorm(100)
ralt = function(mu=0) rnorm(100, mu)
gof_power(pnull, NA, rnull, ralt, c(0, 1), TSextra=TSextra, B=Bsim)

## -----------------------------------------------------------------------------
pnull = function(x, p=c(0,1)) pnorm(x, p[1], ifelse(p[2]>0, p[2], 0.01))
TSextra = list(qnull = function(x, p=c(0,1)) qnorm(x, p[1], ifelse(p[2]>0, p[2], 0.01)))
rnull = function(p=c(0,1)) rnorm(500, p[1], p[2])
ralt = function(mu=0) rnorm(100, mu)
phat = function(x) c(mean(x), sd(x))
gof_power(pnull, NA, rnull, ralt, c(0, 1), phat= phat, 
          TSextra=TSextra, B=Bsim, maxProcessor=2)

## -----------------------------------------------------------------------------
ralt = function(df=1) {
# t distribution truncated at +- 5  
  x=rt(1000, df)
  x=x[abs(x)<5]
  x[1:100]
}  
gof_power(pnull, NA, rnull, ralt, c(2, 50), phat=phat, 
          Range=c(-5,5), TSextra=TSextra, B=Bsim, maxProcessor=2)

## -----------------------------------------------------------------------------
newTScont = function(x, Fx) {
   Fx=sort(Fx)
   n=length(x)
   out = sum(abs( (2*1:n-1)/2/n-Fx ))
   names(out) = "CvM alt"
   out
}

## -----------------------------------------------------------------------------
pnull = function(x) punif(x)
rnull = function() runif(500)
x = rnull()
Rgof::gof_test(x, NA, pnull, rnull, TS=newTScont)

## -----------------------------------------------------------------------------
ralt = function(slope=0) {
  if(slope==0) y=runif(500)
    else y=(slope-1+sqrt((1-slope)^2+4*slope* runif(500)))/2/slope
}

## -----------------------------------------------------------------------------
gof_power(pnull, NA, rnull, ralt, TS=newTScont, param_alt=round(seq(0, 0.5, length=5), 3), Range=c(0,1), B=Bsim)

## -----------------------------------------------------------------------------
vals=1:50/51
pnull = function() (1:50)/50
rnull = function() c(rmultinom(1, 500, rep(1/50,50)))
x = rnull()
gof_test(x, vals, pnull, rnull, TS=newTSdisc)

## -----------------------------------------------------------------------------
ralt = function(slope=0) {
    if(slope==0) p=rep(1/50, 50)
    else p=diff(slope * (0:50/50)^2 + (1 - slope) * 0:50/50)  
  c(rmultinom(1, 500, p))
}
gof_power(pnull, vals, rnull, ralt, TS=newTSdisc, param_alt=round(seq(0, 0.5, length=5), 3), B=Bsim)

## ----eval=FALSE---------------------------------------------------------------
#  pnull=function(x) punif(x)
#  rnull=function() runif(250)
#  pvals=matrix(0,1000,16)
#  for(i in 1:1000) pvals[i, ]=Rgof::gof_test(rnull(), NA, pnull, rnull,B=1000)$p.values

## ----eval=FALSE---------------------------------------------------------------
#  colnames(pvals)=names(Rgof::gof_test(rnull(), NA, pnull, rnull,B=10)$p.values)
#  p1=apply(pvals[, c("W", "ZC", "AD", "ES-s-P" )], 1, min)
#  p2=apply(pvals[, c("KS", "K", "AD", "CvM")], 1, min)

## -----------------------------------------------------------------------------
tmp=readRDS("../inst/extdata/pvaluecdf.rds")
Tests=factor(c(rep("Identical Tests", nrow(tmp)),
        rep("Correlated Selection", nrow(tmp)),
        rep("Best Selection", nrow(tmp)),
        rep("Independent Tests", nrow(tmp))),
        levels=c("Identical Tests",  "Correlated Selection", 
                 "Best Selection", "Independent Tests"),
        ordered = TRUE)
dta=data.frame(x=c(tmp[,1],tmp[,1],tmp[,1],tmp[,1]),
          y=c(tmp[,1],tmp[,3],tmp[,2],1-(1-tmp[,1])^4),
          Tests=Tests)
ggplot2::ggplot(data=dta, ggplot2::aes(x=x,y=y,col=Tests))+
  ggplot2::geom_line(linewidth=1.2)+
  ggplot2::labs(x="p value", y="CDF")+
  ggplot2::scale_color_manual(values=c("blue","red", "Orange", "green"))

## -----------------------------------------------------------------------------
df=3
pnull=function(x) pnorm(x)/(2*pnorm(3)-1)
rnull=function() {x=rt(2000, df);x=x[abs(x)<3];sort(x[1:1000])}
w=function(x) (dnorm(x)/(2*pnorm(3)-1))/(dt(x,df)/(2*pt(3,df)-1))
x=sort(rnull())
plot(x, w(x), type="l", ylim=c(0, 2*max(w(x))))
ralt=function(m=0) {x=rt(2000,df)+m;x=x[abs(x)<3];sort(x[1:1000])}
set.seed(111)
Rgof::gof_power(pnull, NA, rnull, ralt, w=w, param_alt = c(0,0.2), Range=c(-3,3),B=Bsim)


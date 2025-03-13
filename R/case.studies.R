#' This function creates the functions needed to run the various case studies.
#' @param which name of the case study.
#' @param nsample =500, sample size.
#' @return a list of functions

case.studies=function(which, nsample=500) {

  if(which=="uniform.linear.cont") {
    return(list(
      pnull = function(x) punif(x),
      TSextra = list(qnull = function(x) qunif(x)),
      rnull = function() runif(nsample),
      ralt = function(slope) {
         if(slope==0) y=runif(nsample)
         else y=(slope-1+sqrt((1-slope)^2+4*slope* runif(nsample)))/2/slope
         y},
      param_alt=round(seq(0, 0.5, length=25), 3),
      Range=c(0,1)  
    ))
  }  
  if(which=="uniform.linear.disc") {
    return(list(     
      vals=1:50/51,
      pnull = function() (1:50)/50,
      rnull = function() c(rmultinom(1, nsample, rep(1/50,50))),
      ralt = function(slope) {
          if(slope==0) p=rep(1/50, 50)
          else p=diff(slope * (0:50/50)^2 + (1 - slope) * 0:50/50)  
          c(rmultinom(1, nsample, p))},
      param_alt=round(seq(0, 0.5, length=25), 3)
    ))    
  } 
  if(which=="uniform.quadratic.cont") {
    return(list(
      param_alt = round(seq(0, 6, length=25), 3),
      pnull = function(x) punif(x),
     TSextra =  list(qnull = function(x) qunif(x)),
     rnull = function() runif(nsample),
     ralt=function(a) {
       n = nsample
       if(a==0) return(runif(n))
       y=rep(0,n)
       for(i in 1:n) {
         repeat {
         x=runif(1)
         if(runif(1)<(1+a*(x-0.5)^2)/(1+a/4)) {
           y[i]=x
           break
         }
       } 
      }
      y
    },
    param_alt = round(seq(0, 6, length=25), 3),
    Range=c(0,1)    
    ))  
  }
  if(which=="uniform.quadratic.disc") { 
    return(list(
      vals=1:50/51,
      pnull = function() (1:50)/50,
      rnull = function() c(rmultinom(1, nsample, rep(1/50,50))),
      ralt = function(a) {
        if(a==0) p=rep(1/50, 50)
        else p=diff( (24*0:50/50+a+8*a*((0:50/50)-1/2)^3)/(24+2*a) )
      c(rmultinom(1, nsample, p))},
      param_alt = round(seq(0, 6, length=25), 3)
    ))
  }    
  if(which=="uniform.bump.cont") { 
    return(list(
      pnull = function(x) punif(x),
      TSextra = list(qnull = function(x) qunif(x)),
      rnull = function() runif(nsample),
      ralt=function(n) c(runif(nsample-n), rnorm(n, 0.5, 0.05)),
      param_alt=round(seq(0, 150, length=25)),
      Range=c(0,1)    
   ))
  }    
  if(which=="uniform.bump.disc") { 
    return(list(
      vals=1:50/51,
      pnull = function()(1:50)/50,
      rnull = function() c(rmultinom(1, nsample, rep(1/50,50))),
      ralt = function(n) {
        bins = 0:50/50
        p=diff( ( (nsample-n)*punif(bins)+n*pnorm(bins, 0.5, 0.05))/nsample)
        c(rmultinom(1, nsample, p))},
      param_alt=round(seq(0, 150, length=25))
    ))
  }  
  if(which=="uniform.sine.cont") { 
    return(list(
     pnull = function(x) punif(x),
     TSextra = list(qnull = function(x) qunif(x)),
     rnull = function() runif(nsample),
     ralt =  function(a) {
       if(a==0) y = runif(nsample)
       else {
         y = NULL
         repeat { 
           z=runif(nsample)
           I=ifelse(runif(nsample)<(1+a*sin(4*pi*z))/(1+a), TRUE, FALSE)
           y = c(y, z[I]) 
           if(length(y)>nsample) break
         }
       }  
       y=y[1:nsample]
     },
    param_alt=round(seq(0, 0.7, length=25), 3),
    Range=c(0,1)    
   ))
 }
 if(which=="uniform.sine.disc") { 
    return(list(
      vals=1:50/51,
      pnull = function() (1:50)/50,
      rnull = function() c(rmultinom(1, nsample, rep(1/50,50))),
      ralt =  function(a) {
        psin=function(x, a=0) x-a*(cos(4*pi*x)-1)/4/pi
        bins = 0:50/50
        vals = (bins[-1]+bins[-51])/2
        py=diff(psin(bins, a))
        c(rmultinom(1, nsample, py))
      },  
      param_alt=round(seq(0, 0.7, length=25), 3)
    ))
  }   
  if(which=="beta22.betaaa.cont") { 
    return(list(
      pnull = function(x) pbeta(x, 2, 2),
      TSextra = list(qnull = function(x) qbeta(x, 2, 2)),
      rnull = function() rbeta(nsample, 2, 2),
      ralt = function(b=2) rbeta(nsample, b, b),
      param_alt=round(seq(2, 3, length=25), 3),
      Range=c(0,1)    
    ))
  }    
  if(which=="beta22.betaaa.disc") { 
    return(list(
      vals=1:50/51,
      pnull = function() {
        bins=0:50/50
        pbeta(bins[-1], 2, 2)
      }, 
      rnull=function() {
        bins=0:50/50
        p=diff(pbeta(bins, 2, 2))
        c(rmultinom(1, nsample, p))
      },  
      ralt=function(b=2) {
        bins=0:50/50
        p=diff(pbeta(bins, b, b))
        c(rmultinom(1, nsample, p))
      },
      param_alt=round(seq(2, 3, length=25), 3)
    ))
  }    
  if(which=="beta22.beta2a.cont") { 
    return(list(
      pnull = function(x) pbeta(x, 2, 2),
      TSextra = list(qnull = function(x) qbeta(x, 2, 2)),
      rnull = function() rbeta(nsample, 2, 2),
      ralt = function(b=2) rbeta(nsample, 2, b),
      param_alt=round(seq(2, 3, length=25), 3),
      Range=c(0,1)    
    ))
  }    
  if(which=="beta22.beta2a.disc") { 
    return(list(
      vals=1:50/51,
      pnull = function() {
        bins=0:50/50
        pbeta(bins[-1], 2, 2)
      },
      rnull=function() {
        bins=0:50/50
        p=diff(pbeta(bins, 2, 2))
        c(rmultinom(1, nsample, p))
      },
      ralt=function(b=2) {
        bins=0:50/50
        p=diff(pbeta(bins, 2, b))
        c(rmultinom(1, nsample, p))
      },
      param_alt=round(seq(2, 3, length=25), 3)
    ))
  }    
  if(which=="normal.shift.cont") { 
    return(list(
      pnull = function(x) pnorm(x),
      TSextra = list(qnull = function(x) qnorm(x)),
      rnull = function() rnorm(nsample),
      ralt = function(mu=0) rnorm(nsample, mu),
      param_alt=round(seq(0, 0.4, length=25), 3),
      Range=c(-Inf, Inf)
    ))
  }    
  if(which=="normal.shift.disc") { 
    bins=seq(-2.5, 2.5, length=50-1)
    return(list(
      vals=c(-2.8, (bins[-1]+bins[-49])/2, 2.8),
      pnull = function() {
        bins=seq(-2.5, 2.5, length=49)
        c(pnorm(bins), 1)
      },
      rnull=function() {
        bins=seq(-2.5, 2.5, length=49) 
        p=diff(c(0, pnorm(bins), 1))
      c(rmultinom(1, nsample, p))
      },
      ralt=function(mu=0) {
        bins=seq(-2.5, 2.5, length=49) 
        p=diff(c(0, pnorm(bins, mu), 1))
        c(rmultinom(1, nsample, p))
      },
      param_alt=round(seq(0, 0.4, length=25), 3)
    ))
  }    
  if(which=="normal.stretch.cont") { 
    return(list(
      pnull = function(x) pnorm(x),
      TSextra = list(qnull = function(x) qnorm(x)),
      rnull = function() rnorm(nsample),
      ralt = function(s=1) rnorm(nsample, 0, s),
      param_alt=round(seq(1, 1.4, length=25), 3),
      Range=c(-Inf, Inf)
    ))   
  }  
  if(which=="normal.stretch.disc") {  
    bins=seq(-2.5, 2.5, length=50-1) 
    return(list(
    vals=c(-2.8, (bins[-1]+bins[-49])/2, 2.8),
    pnull = function() {
      bins=seq(-2.5, 2.5, length=49)
      c(pnorm(bins), 1)
    },
    rnull=function() {
      bins=seq(-2.5, 2.5, length=49) 
      p=diff(c(0, pnorm(bins), 1))
      c(rmultinom(1, nsample, p))
    },
    ralt=function(s=1) {
      bins=seq(-2.5, 2.5, length=49) 
      p=diff(c(0, pnorm(bins, 0, s), 1))
      c(rmultinom(1, nsample, p))
    },
    param_alt=round(seq(1, 1.4, length=25), 3)
   ))
  }   
  if(which=="normal.t.cont") { 
    return(list(
      pnull = function(x) pnorm(x),
      TSextra = list(qnull = function(x) qnorm(x)), 
      rnull = function() rnorm(nsample),
      ralt = function(df=1) rt(nsample, df),
      param_alt=2*2:26,
      Range=c(-Inf, Inf)
    ))
  }    
  if(which=="normal.t.disc") { 
    bins=seq(-2.5, 2.5, length=49)
    return(list(
      vals=c(-2.8, (bins[-1]+bins[-49])/2, 2.8),
      pnull = function() {
        bins=seq(-2.5, 2.5, length=49)
        c(pnorm(bins), 1)
      },
      rnull=function() {
        bins=seq(-2.5, 2.5, length=49) 
        p=diff(c(0, pnorm(bins), 1))
        c(rmultinom(1, nsample, p))
      },
      ralt=function(df=1) {
        bins=seq(-2.5, 2.5, length=49) 
        p=diff(c(0, pt(bins, df), 1))
        c(rmultinom(1, nsample, p))
      },
      param_alt=2*2:26
    ))
  }    
  if(which=="normal.outlier1.cont") { 
    return(list(  
      pnull = function(x) pnorm(x),
      TSextra = list(qnull = function(x) qnorm(x)),
      rnull = function() rnorm(nsample),
      ralt = function(n=0) c(rnorm(nsample-2*n), runif(2*n, 2, 3)),
      param_alt=0:24,
      Range=c(-Inf, Inf)
   ))
  }   
  if(which=="normal.outlier1.disc") { 
    bins=seq(-3, 3, length=49)
    return(list(
      vals=c(-3.4, (bins[-1]+bins[-49])/2, 3.4),
      pnull = function() {
        bins=c(-Inf, seq(-3, 3, length=49), Inf)
        pnorm(bins[-1])
      },
      rnull=function() {
        bins=c(-Inf, seq(-3, 3, length=49), Inf)
        hist(rnorm(nsample), bins, plot = FALSE)$counts
      },
      ralt=function(n=0) {
        bins=c(-Inf, seq(-3, 3, length=49), Inf)
        hist(c(rnorm(nsample-2*n), runif(2*n, 2, 3)), bins, plot = FALSE)$counts
      },
      param_alt=0:24
    ))
  }    
  if(which=="normal.outlier2.cont") { 
    return(list( 
      pnull = function(x) pnorm(x),
      TSextra = list(qnull = function(x) qnorm(x)),
      rnull = function() rnorm(nsample),
      ralt = function(n=0) c(rnorm(nsample-2*n), runif(n, -3, -2), runif(n, 2, 3)),
      param_alt=0:24
    ))
  }    
  if(which=="normal.outlier2.disc") { 
    bins=seq(-3, 3, length=49)
    return(list(   
      vals=c(-3.4, (bins[-1]+bins[-49])/2, 3.4),
      pnull = function() {
        bins=c(-Inf, seq(-3, 3, length=49), Inf)
        pnorm(bins[-1])
      },
      rnull=function() {
        bins=c(-Inf, seq(-3, 3, length=49), Inf)
        hist(rnorm(nsample), bins, plot = FALSE)$counts
      },
      ralt=function(n=0) {
        bins=c(-Inf, seq(-3, 3, length=49), Inf)
        hist(c(rnorm(nsample-2*n), runif(n, -3, -2), runif(n, 2, 3)), bins, plot = FALSE)$counts
      },
      param_alt=0:24
    ))
  }    
  if(which=="exponential.gamma.cont") { 
    return(list( 
      pnull = function(x) pexp(x, 1),
      TSextra = list(qnull = function(x) qexp(x, 1)),
      rnull = function() rexp(nsample, 1),
      ralt = function(b=1) rgamma(nsample, 1, b),
      param_alt=round(seq(1, 1.3, length=25), 3),
      Range=c(0, Inf)
    ))
  }    
  if(which=="exponential.gamma.disc") { 
    bins = seq(0, 4, length=50)
    return(list( 
      vals = c((bins[-1]+bins[-50])/2, 5),
      pnull = function() {
        bins = seq(0, 4, length=50)
        pexp(c(bins[-1], Inf), 1)
      },  
      rnull=function() {
        bins = seq(0, 4, length=50)
        p=diff(pexp(c(bins, Inf), 1))  
        c(rmultinom(1, nsample, p))
      },
      ralt=function(b=1) {
        bins = seq(0, 4, length=50)
        p=diff(pgamma(c(bins, Inf),  1, b))  
        c(rmultinom(1, nsample, p))
      },
      param_alt=round(seq(1, 1.3, length=25), 3)
    ))
   }   
  if(which=="exponential.weibull.cont") { 
    return(list(
      pnull = function(x) pexp(x, 1),
      TSextra = list(qnull = function(x) qexp(x, 1)),
      rnull = function() rexp(nsample, 1),
      ralt = function(b=1) rweibull(nsample, 1, b),
      param_alt=round(seq(1, 1.5, length=25), 3),
      Range=c(0,Inf)
    ))
  }    
  if(which=="exponential.weibull.disc") { 
    bins = seq(0, 4, length=51)
    return(list(
      vals = (bins[-1]+bins[-51])/2,
      pnull = function() {
        bins = seq(0, 4, length=50)
        pexp(c(bins[-1], Inf), 1)
      },    
      rnull=function() {
        bins = seq(0, 4, length=50)
        p=diff(pexp(c(bins, Inf), 1))  
        c(rmultinom(1, nsample, p))
      },
      ralt=function(b=1) {
        bins = seq(0, 4, length=50)
        p=diff(pweibull(c(bins, Inf),  1, b))  
        c(rmultinom(1, nsample, p))
      },
      param_alt=round(seq(1, 1.5, length=25), 3)
    ))
  }    
  if(which=="exponential.bump.cont") { 
    return(list(
      pnull = function(x) pexp(x, 1),
      TSextra = list(qnull = function(x) qexp(x, 1)),
      rnull = function() rexp(nsample, 1),
      ralt = function(n=0) c(rexp(nsample-n, 1), rnorm(n, 0.5, 0.05)),
      param_alt=3*0:24,
      Range=c(0, Inf)
    ))
  }    
  if(which=="exponential.bump.disc") { 
    bins = seq(0, 4, length=51)
    return(list(
      vals = (bins[-1]+bins[-51])/2,
      pnull = function() {
        bins = seq(0, 4, length=50)
        pexp(c(bins[-1], Inf), 1)
      },  
      rnull=function() {
        bins = seq(0, 4, length=50)
        p=diff(pexp(c(bins, Inf), 1))  
        c(rmultinom(1, nsample, p))
      },
      ralt=function(n=0) {
        bins = c(seq(0, 4, length=50), Inf)
        p=diff( ((nsample-n)*pexp(bins, 1)+n*pnorm(bins, 0.5, 0.05))/nsample )
        c(rmultinom(1, nsample, p))
      },
      param_alt=3*0:24
    ))
  }    
  if(which=="trunc.exponential.linear.cont") { 
    return(list(
      pnull = function(x) (1-exp(-x))/(1-exp(-1)),
      TSextra = list(qnull = function(x) -log(1 - x*(1-exp(-1)))),
      rnull = function() {
          x <- NULL
          repeat {
              x = c(x, rexp(nsample, 1))
              x = x[x < 1]
              if (length(x) > nsample) 
                  break
          }
          x[1:nsample]
      },
      ralt = function(slope) {
        if(slope==0) y=runif(nsample)
        else y=(slope-1+sqrt((1-slope)^2+4*slope* runif(nsample)))/2/slope
      },
      param_alt=round(seq(-0.5, -1, length=25), 3),
      Range=c(0,1)
    ))
  }    
  if(which=="trunc.exponential.linear.disc") { 
    bins=0:50/50
    return(list(
      vals=(bins[-1]+bins[-51])/2,
      pnull = function() {
        bins=1:50/50
        (1-exp(-bins))/(1-exp(-1))
      },
      rnull=function() {
        bins=0:50/50  
        p=diff((1-exp(-bins))/(1-exp(-1)))
        c(rmultinom(1, 1000, p))
      },
      ralt = function(slope) {
          if(slope==0) p=rep(1/50, 50)
          else p=diff(slope * (0:50/50)^2 + (1 - slope) * 0:50/50)  
        c(rmultinom(1, 1000, p))
      },
      param_alt=round(seq(-0.5, -1, length=25), 3)
    ))
  }    
  if(which=="normal.t.est.cont") { 
    return(list(
      pnull = function(x, p=c(0,1)) pnorm(x, p[1], p[2]),
      TSextra=list(qnull = function(x, p=c(0,1)) qnorm(x, p[1], p[2])),
      rnull = function(p=c(0,1)) rnorm(nsample, p[1], p[2]),
      ralt = function(df=1) {
        x=NULL
        repeat {
          x=c(x, rt(nsample, df))
          x=x[abs(x)<10]
          if(length(x)>=nsample) break
        }  
        x[1:nsample]
      },
      phat = function(x) c(mean(x), sd(x)),
      param_alt=1:25,
      Range=c(-Inf, Inf)
    ))
  }    
  if(which=="normal.t.est.disc") { 
    bins = c(-2.5, seq(-2.5, 2.5, length=49), 2.5)
    return(list(
      vals = (bins[-1]+bins[-51])/2,
      pnull = function(p=c(0,1)) {
        bins = c(seq(-2.5, 2.5, length=49), Inf)
        pnorm(bins, p[1], abs(p[2]))
      },
      rnull=function(p=c(0,1)) {
        bins = c(-Inf, seq(-2.5, 2.5, length=49), Inf)
        p=diff(pnorm(bins, p[1], abs(p[2])))
        c(rmultinom(1, nsample, p))
      },
      ralt=function(df=1) {
        bins = c(-Inf, seq(-2.5, 2.5, length=49), Inf)
        p=diff(pt(bins, df))
        c(rmultinom(1, nsample, p))
      },
      phat = function(x) {
        bins = c(-2.5, seq(-2.5, 2.5, length=49), 2.5)
        vals = (bins[-1]+bins[-51])/2
        c(mean(rep(vals, x)), sd(rep(vals, x)))
      },
      param_alt=1:25
    ))
  }    
  if(which=="exponential.weibull.est.cont") { 
    return(list(
      pnull = function(x, p=1) pexp(x, p),
      TSextra=list(qnull = function(x, p=1) qexp(x, p)),
      rnull = function(p=1) rexp(nsample, p),
      ralt = function(b=1) rweibull(nsample, b, 1),
      phat = function(x) 1/mean(x),
      param_alt = round(seq(1, 1.4, length=25), 3),
      Range=c(0, Inf)
    ))
  }    
  if(which=="exponential.weibull.est.disc") { 
    bins=seq(0, 4, length=51)
    return(list(
      vals=(bins[-1]+bins[-51])/2,
      pnull = function(a=1) {
        bins=seq(0, 4, length=51)
        pexp(c(bins[2:50],Inf), a)
      },
      rnull=function(a=1) {
        bins=seq(0, 4, length=51)
        p = diff(pexp(c(0, bins[2:50],Inf), a))
        c(rmultinom(1, nsample, p))
      },
      ralt=function(b=1) {
        bins=seq(0, 4, length=51)
        p = diff(pweibull(c(0, bins[2:50], Inf), b, 1))
        c(rmultinom(1, nsample, p))
      },  
      phat = function(x) {
        bins=seq(0, 4, length=51)
        vals=(bins[-1]+bins[-51])/2
        1/mean(rep(vals, x))  
      },
      param_alt = round(seq(1, 1.4, length=25), 3)
    ))
  }    
  if(which=="trunc.exponential.linear.est.cont") { 
    return(list(
      pnull = function(x, p=1) (1-exp(-x*p))/(1-exp(-p)),
      TSextra=list(qnull = function(x, p=1) -log(1 - x*(1-exp(-p)))/p),
      rnull = function(p=1) {
          x = NULL
          repeat {
            x = c(x, rexp(1000, p))
            x = x[x < 1]
            if (length(x) > 1000) 
              break
          }
          x[1:1000] 
      },
      ralt = function(slope) {
        if(slope==0) y=runif(1000)
        else y=(slope-1+sqrt((1-slope)^2+4*slope* runif(1000)))/2/slope
      },
      phat = function(x) {
          n = length(x)
          s = sum(x)
          p = n/s
          repeat {
          o = p
            t = exp(-o)
            l1 =  n/o - s - n*t/(1-t)
            l2 = (-n/o^2 + n * t/(1-t)^2)
            p = o - l1/l2
            if(p<0) return(0.001)
            if (abs(p - o) < 0.01) 
                break
          } 
          p
      },
      param_alt=round(seq(-0.5, -1, length=25), 3),
      Range=c(0, 1)
    ))
  }    
  if(which=="trunc.exponential.linear.est.disc") { 
    bins = 0:100/100
    return(list(
      vals = (bins[-1]+bins[-101])/2,
      pnull = function(a=1) {
        bins=1:100/100
        (1-exp(-bins*a))/(1-exp(-a))
      },
      rnull=function(a=1) {
        bins=0:100/100
        p=diff(1-exp(-bins*a))
        c(rmultinom(1, 1000, p))
      },
      ralt = function(slope) {
        if(slope==0) p=rep(1/100, 100)
        else p=diff(slope * (0:100/100)^2 + (1 - slope) * 0:100/100)  
        c(rmultinom(1, 1000, p))
      },
      phat = function(x) {
        bins=0:100/100
        vals=(bins[-1]+bins[-101])/2
        n = sum(x)
        s = sum(vals*x)
        p = n/s
        repeat {
            o = p
            t = exp(-o)
            l1 =  n/o - s - n*t/(1-t)
            l2 = (-n/o^2 + n * t/(1-t)^2)
            p = o - l1/l2
            if(p<0) return(0.01)
            if (abs(p - o) < 0.01) 
                break
        }
        p
     },
     param_alt=round(seq(-0.5, -1, length=25), 3)
   ))
  }
  if(which=="exponential.gamma.est.cont") { 
    return(list(
      pnull = function(x, p) pexp(x, p),
      TSextra=list(qnull = function(x, p) qexp(x, p)),
      rnull = function(p) rexp(nsample, p),
      ralt = function(b=1) rgamma(nsample, b, 1),
      phat = function(x) 1/mean(x),
      param_alt=round(seq(1, 1.3, length=25), 3),
      Range=c(0, Inf)
    ))
  }    
  if(which=="exponential.gamma.est.disc") { 
    bins = seq(0, 4, length=51)
    return(list(
      vals = (bins[-1]+bins[-51])/2,
      pnull = function(b=1) {
        bins = seq(0, 4, length=50)
        pexp(c(bins[-1], Inf), b)
      },  
      rnull=function(b=1) {
        bins = seq(0, 4, length=50)
        p=diff(pexp(c(bins, Inf), b))  
        c(rmultinom(1, nsample, p))
      },
      ralt=function(b=1) {
        bins = seq(0, 4, length=50)
        p=diff(pgamma(c(bins, Inf), b, 1))
        c(rmultinom(1, nsample, p))
      },
      phat=function(x) {
        bins = seq(0, 4, length=51)
        vals = (bins[-1]+bins[-51])/2
        1/mean(rep(vals, x))
      },
      param_alt=round(seq(1, 1.3, length=25), 3),
      Range=c(0, Inf)
    ))
  }    
  if(which=="normal.cauchy.est.cont") { 
    return(list(
      pnull = function(x, s=1) pnorm(x, 0, s),
      TSextra=list(qnull = function(x, s=1) qnorm(x, 0, s)),
      rnull = function(s=1) rnorm(nsample, 0, s),
      ralt = function(s) s*tan( (2*runif(nsample)-1)*atan(5/s)),
      phat = function(x) {
        ll=function(s, x) -sum(log(dnorm(x,0,s))-log(2*pnorm(5,0,s)-1))
        optimize(ll, c(0,20), x=x)$minimum
      },
      param_alt=round(seq(1, 2, length=25), 3),
      Range=c(-Inf, Inf)
   ))
 }     
  if(which=="normal.cauchy.est.disc") { 
    bins = seq(-5, 5, 0.25)
    return(list(
      vals = (bins[-1]+bins[-41])/2,
      pnull = function(s=1) {
        bins = seq(-5, 5, 0.25)
        pnorm(bins[-1], 0, s)/(2*pnorm(5, 0, s)-1)
      },   
      rnull=function(s=1) {
        bins=seq(-5, 5, 0.25)
        p=diff(pnorm(bins, 0, s))
        c(rmultinom(1, nsample, p))
      },
      ralt = function(s=1) {
        bins=seq(-5, 5, 0.25)
        p=diff(pcauchy(bins, 0, s))
        c(rmultinom(1, nsample, p))
      },   
      phat = function(x) {
        ll=function(s, x) {
          bins = seq(-5, 5, 0.25)
          vals = (bins[-1]+bins[-41])/2       
              -sum(x*log(dnorm(vals,0,s)))+sum(x)*log(2*pnorm(5,0,s)-1)
        }    
        optimize(ll, c(0, 4), x=x)$minimum
      },
      param_alt=round(seq(1, 2, length=25), 3)
    ))
  }  
}

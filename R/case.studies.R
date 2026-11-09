#' Create case studies
#' 
#' This function creates the functions needed to run the various case studies.
#'
#' @param which name of the case study, or its number in the list.
#' @param which_diff ="mean", name of parameter in which data sets differ, 
#'                    or its number in the list.
#' @param data_type ="cont", is data continuous or discrete?
#' @param n sample sizes.
#' @param nbins =50, number of bins for discrete data
#' @param funs_list list of alternative functions
#' @param ReturnCaseNames =FALSE, should list of case studies be returned?
#' @return a list of functions to create case studies
#' @export 
case.studies=function(which, which_diff ="mean", data_type="cont",
                n, nbins=50, funs_list, ReturnCaseNames =FALSE) {
  # check if old version is assumed
  if(length(which_diff)==1 &&
     is.numeric(which_diff) &&
     !(which_diff %in% 1:6)) {
    stop("the arguments list of the function case.studies has changed.
          Please check the help file or the vignette",
         call.=FALSE)
  }

  list.of.studies=c(
    "uniform", "linear", "quadratic", "uniform-betabump", "sine", "betaaa",
    "beta2a", "triangular","uniformmixture", "betamixture", "normal-pure",
    "normal-mix1",  "normal-mix2", "normal-mix3", "normal-mix4", "gamma1", 
    "gamma2", "gammamixture",  "exponential-mixture", "noncentral",
    "uniform-linear", "uniform-quadratic", "uniform-bump",
    "uniform-sine",  "beta22.betaaa", "beta22.beta2a", 
    "uniform.uniformmixture", "uniform.betamixture",
    "uniform.triangular", "double.exponential", "normal.onestretch",
    "normal.t", "normal.outlier1", 
    "normal.outlier2", "normal.normalmixture",
    "exponential.gamma", "exponential.weibull", 
    "exponential.bump", "chisquare.noncentral",  
    "gamma.gammamixture", "normal.t.est", "exponential.weibull.est",
    "trunc.exponential.linear.est", "exponential.gamma.est", 
    "normal.cauchy.est",
    "betaa1.betaab.est", "betaa1.betaaa.est",
    "trans-beta.est", "double.exponential.est",
    "pareto.est", "betamixture1.est", 
    "betamixture2.est", "betamixture3.est", "normalmixture.est",
    "laplace1.est", "laplace2.est","truncexp1.est",
    "truncexp2.est", "gammamixture1.est", "gammamixture2.est"
  )
  if(ReturnCaseNames) return(list.of.studies)

  if(is.numeric(which[1])) which=list.of.studies[which]
  alltypes=c("mean", "variance", "skewness", "kurtosis", "mixed")
  if(is.numeric(which_diff[1]))
    which_diff=alltypes[which_diff]
  if(missing(n)) {
    n=Rgof::case_studies_sample_sizes[[data_type]]
    if(which %in% list.of.studies[1:20]) 
         n=n[[1]][which, which_diff]
    else n=n[[2]][which]
  }
  if(which %in% list.of.studies[1:20]) {
    if(!(which_diff%in%alltypes[1:4]))
       stop("for case studies #1-20 which_diff must be mean, variance, skewness, or kurtosis, or 1-4",
         call.=FALSE)
    if(missing(funs_list)) {
      funs = Rgof::funs_list[[data_type]][[which]][[which_diff]]
    } else {
      funs = funs_list
    }
    apar=0
    ralt = function(a=apar[1]) funs$random(n)
    palt = function(x, a=apar[1]) funs$cdf(x)
  }
  
  phat=function(x) -99
  # quantile function by vectorized bisection: f is called about 50 times with
  # the whole vector p, instead of uniroot once for each element of p
  finv <- function(p, f, Ra) {
    lo <- rep(Ra[1], length(p))
    hi <- rep(Ra[2], length(p))
    for(k in 1:50) {
      mid <- (lo+hi)/2
      below <- f(mid) < p
      lo[below] <- mid[below]
      hi[!below] <- mid[!below]
    }
    y <- (lo+hi)/2
    y[p <= 0] <- -Inf
    y[p >= 1] <- Inf
    if (f(0) == 0) y[p == 0] <- 0
    if (f(1) == 1) y[p == 1] <- 1
    y
  }
  # quantile function of a null distribution without parameters, from a table:
  # F is evaluated once (on first use) at m equally spaced points of the range Ra
  # and at points close to its ends, then the inverse is interpolated with a
  # monotone cubic spline. Much faster than finv; the error is below 1e-6, except
  # for probabilities below 1e-6 or above 1-1e-6 (below 1e-4 there).
  qtab <- function(f, Ra, m=20001) {
    Q <- NULL
    lo <- NULL
    hi <- NULL
    function(p) {
      if(is.null(Q)) {
        d <- diff(Ra)*10^seq(-12, -4, length.out=200)
        x <- sort(unique(c(seq(Ra[1], Ra[2], length.out=m), Ra[1]+d, Ra[2]-d)))
        Fx <- f(x)
        k <- c(TRUE, diff(Fx) > 0)
        Q <<- stats::splinefun(Fx[k], x[k], method="monoH.FC")
        lo <<- c(min(Fx), if(Fx[1] == 0) Ra[1] else -Inf)
        hi <<- c(max(Fx), if(Fx[length(x)] == 1) Ra[2] else Inf)
      }
      y <- Q(pmin(pmax(p, lo[1]), hi[1]))
      y[p <= 0] <- lo[2]
      y[p >= 1] <- hi[2]
      y
    }
  }
# sizes of mixture components that add up to n exactly
  msizes=function(n, p) {
    k=floor(n*p)
    k[length(k)]=n-sum(k[-length(k)])
    k
  }
  collect_case <- function(env) {
    keep <- c("pnull", "dnull", "rnull", "palt", "ralt", "TSextra",
              "phat", "phatD", "pnullD", "rnullD", "bins", "Ra", "apar")
    present <- vapply(keep, exists, logical(1), envir=env, inherits=FALSE)
    mget(keep[present], envir=env, inherits=FALSE)
  }
  
  study <- switch(
    which,
    # 1
    "uniform" = (function() {
      pnull=function(x) punif(x)
      dnull=function(x) dunif(x)
      rnull=function() runif(n)
      TSextra=list(qnull=function(x) qunif(x))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 2
    "linear" = (function() {
      pnull=function(x) {
        ifelse(x<0, 0, ifelse(x>1, 1, 0.3*x^2+0.7*x))
      }
      rnull=function() 10/6*(-0.7+sqrt(0.49+1.2*runif(n)))
      dnull=function(x) 
        ifelse(x<0, 0, ifelse(x>1, 0, 0.6*x+0.7))
      TSextra =  list(qnull = function(x) 
        (-7/6+sqrt((7/6)^2+10/3*x)))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 3
    "quadratic" = (function() {
      pnull=function(x) {
        a=-0.5;
        ifelse(x<0, 0, ifelse(x>1, 1, a*(x-1/2)^3+(1-a/4)*x+a/8))
        }
      rnull=function() {
        y=rep(0, n)
        for(i in 1:n) {
           theta <- acos(4 * (1 - 2*runif(1)) / (3*sqrt(3))) / 3
           y[i]=(1/2 + sqrt(3) * cos(theta - 2*pi*(0:2)/3))[2]
        }
        y
      }  
      dnull=function(x) {
        a=-0.5;
        ifelse(x<0, 0, ifelse(x>1, 0, 3*a*(x-1/2)^2+1-a/4))
        }
      TSextra =  list(qnull = qtab(pnull, c(0, 1)))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 4
    "uniform-betabump" = (function() {
      pnull=function(x) 0.825*punif(x)+0.175*pbeta(x, 25, 25)
      rnull=function() c(runif(floor(n*0.825)), 
                          rbeta(ceiling(n*0.175), 25, 25))
      dnull=function(x) 0.825*dunif(x)+0.175*dbeta(x, 25, 25)
      TSextra =  list(qnull = qtab(pnull, c(0, 1)))
      bins=seq(0, 1, length=nbins+1)    
      collect_case(environment())
    })(),
    # 5
    "sine" = (function() {
      pnull=function(x) 
        ifelse(x<0, 0, ifelse(x>1, 1, (1-cos(pi*x))/2))
      rnull=function() acos(2*runif(n)-1)/pi   
      dnull=function(x) 
        ifelse(x<0, 0, ifelse(x>1, 0, pi*sin(pi*x)/2))
      TSextra=list(qnull = qtab(pnull, c(0, 1)))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 6
    "betaaa" = (function() {
      pnull=function(x) pbeta(x, 3, 3)
      rnull=function() rbeta(n, 3, 3)
      dnull=function(x) dbeta(x, 3, 3)
      TSextra=list(qnull = function(x) qbeta(x, 3, 3))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 7
    "beta2a" = (function() {
      pnull=function(x) pbeta(x, 2, 2.5)    
      rnull=function() rbeta(n, 2, 2.5)
      dnull=function(x) dbeta(x, 2, 2.5)
      TSextra=list(qnull = function(x) qbeta(x, 2, 2.5))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 8
    "triangular" = (function() {
       pnull=function(x) {
          y=0.65*punif(x)+0.35*ifelse(x<0.5, 2*x^2, 1-2*(1-x)^2)
          ifelse(x<0, 0, ifelse(x>1, 1, y))
         }
       rnull=function() {
         k=msizes(n, c(0.35, 0.65))
         y=runif(k[1])
         y=ifelse(y<1/2, sqrt(y/2), 1-sqrt((1-y)/2))
         c(y, runif(k[2]))
      }    
      dnull=function(x) {
        y=0.65*dunif(x)+0.35*ifelse(x<0.5, 4*x, 4-4*x)
        ifelse(x<0, 0, ifelse(x>1, 0,y)) 
        }
      TSextra=list(qnull=function(x) ifelse(x<1/2, (-0.65+sqrt(0.65^2+8*0.35*x))/4/0.35, 
                     (2.05-sqrt(2.05^2-4*0.7*(0.35+x)))/1.4))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 9
    "uniformmixture" = (function() {
      pnull=function(x) (punif(x, 0, 1/3)+2*punif(x,1/3,2/3)+3*punif(x, 2/3, 1))/6
      rnull=function() {
        k=msizes(n, c(1/6, 1/3, 1/2))
        c(runif(k[1], 0, 1/3), runif(k[2],1/3,2/3), runif(k[3], 2/3, 1))
      }
      dnull=function(x) (dunif(x, 0, 1/3)+2*dunif(x,1/3,2/3)+3*dunif(x, 2/3, 1))/6
      TSextra=list(
          qnull=function(x) {
             y=2*x
             y[x>1/6]=x[x>1/6]+1/6
             y[x>1/2]=(2*x[x>1/2]+1)/3
             y[x<0]=0
             y[x>1]=1
             y
          }
        )
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 10
    "betamixture" = (function() {
      pnull=function(x) (punif(x)+pbeta(x, 2, 2))/2
      rnull=function() {k=msizes(n, c(1/2, 1/2)); c(runif(k[1]), rbeta(k[2], 2, 2))}
      dnull=function(x) (dunif(x)+dbeta(x, 2, 2))/2
      TSextra=list(qnull =qtab(pnull, c(0, 1)))
      bins=seq(0, 1, length=nbins+1)    
      collect_case(environment())
    })(),
    # 11
    "normal-pure" = (function() {
      pnull=function(x) pnorm(x)
      rnull=function() rnorm(n)
      dnull=function(x) dnorm(x)
      TSextra=list(qnull=function(x) qnorm(x))
      bins=seq(-3, 3, length=nbins+1)  
      bins[1]=-Inf;bins[nbins+1]=Inf
      collect_case(environment())
    })(),
    # 12
    "normal-mix1" = (function() {
      pnull=function(x) (pnorm(x)+pnorm(x, 0, 3))/2
      rnull=function() c(rnorm(floor(n/2)), rnorm(ceiling(n/2), 0, 3))
      dnull=function(x) (dnorm(x)+dnorm(x, 0, 3))/2
      TSextra=list(qnull=qtab(pnull, c(-20, 20)))
      bins=seq(-12, 12, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 13
    "normal-mix2" = (function() {
      pnull=function(x) (pnorm(x)+pnorm(x, 1))/2
      rnull=function() c(rnorm(floor(n/2)), rnorm(ceiling(n/2), 1))
      dnull=function(x) (dnorm(x)+dnorm(x, 1))/2
      TSextra=list(qnull=qtab(pnull, c(-20, 20)))
      bins=seq(-12, 12, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf  
      collect_case(environment())
    })(),
    # 14
    "normal-mix3" = (function() {
      pnull=function(x) pnorm(x, -3)/4+pnorm(x)/2+pnorm(x, 5)/4
      rnull=function() {k=msizes(n, c(1/4, 1/2, 1/4)); c(rnorm(k[1], -3), rnorm(k[2]), rnorm(k[3],5))}
      dnull=function(x) dnorm(x, -3)/4+dnorm(x)/2+dnorm(x, 5)/4
      TSextra=list(qnull=qtab(pnull, c(-20, 20)))
      bins=seq(-3, 4, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 15
    "normal-mix4" = (function() {
      pnull=function(x) (2*pnorm(x)+pnorm(x, 0, 3)+pnorm(x, 0, 5))/4
      rnull=function() {k=msizes(n, c(1/2, 1/4, 1/4)); c(rnorm(k[1]), rnorm(k[2], 0, 3), rnorm(k[3], 0, 5))}
      dnull=function(x) (2*dnorm(x)+dnorm(x, 0, 3)+dnorm(x, 0, 5))/4
      TSextra=list(qnull=qtab(pnull, c(-50, 50)))
      bins=seq(-25, 20, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 16
    "gamma1" = (function() {
      pnull=function(x) 0.2*pexp(x, 1)+0.8*pgamma(x, 2, 0.5)    
      rnull=function() {k=msizes(n, c(0.2, 0.8)); c(rexp(k[1], 1), rgamma(k[2], 2, 0.5))}
      dnull=function(x) 0.2*dexp(x, 1)+0.8*dgamma(x, 2, 0.5)
      TSextra=list(qnull=qtab(pnull, c(0, 100)))
      bins=seq(0, 25, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 17
    "gamma2" = (function() {
      pnull=function(x) 0.2*pgamma(x, 2, 2)+0.8*pgamma(x, 3, 5)    
      rnull=function() {k=msizes(n, c(0.2, 0.8)); c(rgamma(k[1], 2, 2), rgamma(k[2], 3, 5))}
      dnull=function(x) 0.2*dgamma(x, 2, 2)+0.8*dgamma(x, 3, 5)
      TSextra=list(qnull=qtab(pnull, c(0, 40)))
      bins=seq(0, 5, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 18
    "gammamixture" = (function() {
      pnull=function(x) (pgamma(x, 2, 1)+pgamma(x, 10,1))/2
      rnull=function() {k=msizes(n, c(1/2, 1/2)); c(rgamma(k[1], 2, 1), rgamma(k[2], 10,1))}
      dnull=function(x) (dgamma(x, 2, 1)+dgamma(x, 10,1))/2
      TSextra=list(qnull=qtab(pnull, c(0, 100)))
      bins=seq(0, 30, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 19
    "exponential-mixture" = (function() {
      pnull=function(x) (pexp(x, 1)+pexp(x, 3))/2
      rnull=function() c(rexp(floor(n/2), 1), 
              rexp(ceiling(n/2), 3))
      dnull=function(x) (dexp(x, 1)+dexp(x, 3))/2
      TSextra=list(qnull=qtab(pnull, c(0, 30)))
      bins=seq(0, 8, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 20
    "noncentral" = (function() {
      pnull=function(x) 
        ifelse(x< 25, pchisq(x, 5, 0.875)/pchisq(25, 5, 0.875), 1)
      dnull=function(x) 
        ifelse(x< 25, dchisq(x, 5, 0.875)/pchisq(25, 5, 0.875), 0)
      rnull=function() {
        x=rchisq(1.1*n, 5, 0.875)
        x=x[x<25]
        x[1:n]
      }
      TSextra=list(
        qnull=qtab(pnull, c(0, 25))
      )
      bins=seq(0, 25, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 21
    "uniform-linear" = (function() {
      apar=c(0.3, 0)
      rnull=function() runif(n)
      dnull=function(x) dunif(x)
      pnull=function(x) punif(x)
      TSextra=list(qnull=function(x) qunif(x))
      ralt=function(a=apar[1]) {
        if(a==0) return(runif(n))
        (-(1-a)+sqrt((1-a)^2+4*a*runif(n)))/2/a
        }
      palt=function(x, a=apar[1]) {
        if(a==0) return(punif(x))
        ifelse(x<0, 0,
               ifelse(x>1, 1,
                      a*x^2+(1-a)*x))
        }
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 22
    "uniform-quadratic" = (function() {
      apar=ifelse(data_type=="cont", -1.2, -1)
      apar=c(apar, 0)
      rnull=function() runif(n)
      dnull=function(x) dunif(x)
      pnull=function(x) punif(x)
      ralt=function(a=apar[1]) {
        if(a==0) return(runif(n))
        x=rep(0,n)
        for(i in 1:n) {
          repeat {
            z=runif(1)
            if(max(1-a/4, 1+a/2)*runif(1)<3*a*(z-1/2)^2+1-a/4) { # divide by max of density
              x[i]=z
              break
            }
          } 
        }
        x
      }   
      palt=function(x, a=apar[1]) {
        if(a==0) return(punif(x))
        ifelse(x<0, 0,
               ifelse(x>1, 1,
                      a*(x-0.5)^3+x-a*x/4+a/8))
        }
      TSextra=list(qnull=function(x) qunif(x))
      bins=seq(0, 1, length=nbins+1)
      collect_case(environment())
    })(),
    # 23
    "uniform-bump" = (function() {
      apar=c(0.125, 0)
      rnull=function() runif(n)
      dnull=function(x) dunif(x)
      pnull=function(x) punif(x)
      palt=function(x, a=apar[1]) (1-a)*punif(x)+a*pnorm(x, 0.5, 0.05)
      ralt=function(a=apar[1]) c(runif(floor(n*(1-a))), 
          rnorm(ceiling(n*a), 0.5, 0.05))
      TSextra=list(qnull=function(x) qunif(x))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 24
    "uniform-sine" = (function() {
      apar=ifelse(data_type=="cont", 0.5, 0.4)
      apar=c(apar, 0)
      rnull=function() runif(n)
      dnull=function(x) dunif(x)
      pnull=function(x) punif(x)
      palt=function(x, a=apar[1]) {
        if(a==0) return(punif(x))
        ifelse(x<0, 0,
               ifelse(x>1, 1,
                      x-a*cos(4*pi*x)/4/pi+a/4/pi))
        }
      ralt=function(a=apar[1]) {
        if(a==0) return(runif(n))
        y = NULL
        repeat { 
          z=runif(n)
          I=ifelse((1+abs(a))*runif(n)<1+a*sin(4*pi*z), TRUE, FALSE) # divide by max of density
          y = c(y, z[I]) 
          if(length(y)>n) break
        }
        y[1:n]
      }
      TSextra=list(qnull=function(x) qunif(x))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 25
    "beta22.betaaa" = (function() {
      apar=ifelse(data_type=="cont", 2.7, 2.85)
      apar=c(apar, 2)
      rnull=function() rbeta(n, 2, 2)
      dnull=function(x) dbeta(x, 2, 2)
      pnull=function(x) pbeta(x, 2, 2)
      palt=function(x, a=apar[1]) pbeta(x, a, a)
      ralt=function(a=apar[1]) rbeta(n, a, a)
      TSextra=list(qnull=function(x) qbeta(x, 2, 2))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 26
    "beta22.beta2a" = (function() {
      apar=c(2.35, 2)
      rnull=function() rbeta(n, 2, 2)
      dnull=function(x) dbeta(x, 2, 2)
      pnull=function(x) pbeta(x, 2, 2)
      palt=function(x, a=apar[1]) pbeta(x, 2, a)
      ralt=function(a=apar[1]) rbeta(n, 2, a)
      TSextra=list(qnull=function(x) qbeta(x, 2, 2))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 27
    "uniform.uniformmixture" = (function() {
      apar=ifelse(data_type=="cont", 0.57, 0.55)
      apar=c(apar, 0.5)
      rnull=function() runif(n)
      dnull=function(x) dunif(x)
      pnull=function(x) punif(x)
      palt=function(x, a=apar[1]) 
        a*punif(x, 0, 0.5)+(1-a)*punif(x, 0.5, 1)
      ralt=function(a=apar[1]) 
        ifelse(runif(n)<a, runif(n, 0, 0.5), runif(n, 0.5, 1))
      TSextra=list(qnull=function(x) qunif(x))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 28
    "uniform.betamixture" = (function() {
      apar=c(0.5, 0)
      rnull=function() runif(n)
      dnull=function(x) dunif(x)
      pnull=function(x) punif(x)
      palt=function(x, a=apar[1]) 
        (1-a)*punif(x)+a*pbeta(x, 2, 2)
      ralt=function(a=apar[1]) 
        c(runif(floor((1-a)*n)), rbeta(ceiling(a*n), 2, 2))
      TSextra=list(qnull=function(x) qunif(x))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 29
    "uniform.triangular" = (function() {
      apar=c(0.44, 0)
      rnull=function() runif(n)
      dnull=function(x) dunif(x)
      pnull=function(x) punif(x)
      palt=function(x, a=apar[1]) {
        if(a==0) return(punif(x))
        ifelse(x<0, 0,
               ifelse(x>1, 1,
                      (1-a)*x+a*ifelse(x<1/2, 2*x^2, 4*x-2*x^2-1)))
        }
      ralt=function(a=apar[1]) {
        y=runif(ceiling(a*n))
        y=ifelse(y<1/2, sqrt(y/2), 1-sqrt((1-y)/2))
        c(y, runif(floor((1-a)*n)))
      }    
      TSextra=list(qnull=function(x) qunif(x))
      bins=seq(0, 1, length=nbins+1)  
      collect_case(environment())
    })(),
    # 30
    "double.exponential" = (function() {
      apar=ifelse(data_type=="cont", 1.2, 1.3)
      apar=c(apar, 1)
      rnull=function() ifelse(runif(n)<0.5,-rexp(n, 1), rexp(n, 1))
      dnull=function(x) dexp(abs(x), 1)/2
      pnull=function(x) ifelse(x<0,(1-pexp(-x, 1))/2,1/2+pexp(x, 1)/2)
      palt=function(x, a=apar[1]) ifelse(x<0,(1-pexp(-x, a))/2,1/2+pexp(x, a)/2)
      ralt=function(a=apar[1]) ifelse(runif(n)<0.5,-rexp(n, a), rexp(n, a))
      TSextra=list(qnull=function(x) ifelse(x<1/2, log(2*x), -log(2*(1-x))))
      bins=seq(-10, 10, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf
      collect_case(environment())
    })(),
    # 31
    "normal.onestretch" = (function() {
      apar=ifelse(data_type=="cont", 1.25, 1.35) 
      apar=c(apar, 1)
      rnull=function() rnorm(n)
      dnull=function(x) dnorm(x)
      pnull=function(x) pnorm(x)
      palt=function(x, a=apar[1]) ifelse(x<0, pnorm(x), pnorm(x, 0, a))
      ralt=function(a=apar[1]) {k=msizes(n, c(1/2, 1/2)); c(-abs(rnorm(k[1])), abs(rnorm(k[2], 0,a)))}
      TSextra=list(qnull=function(x) qnorm(x))
      bins=seq(-4, 4, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 32
    "normal.t" = (function() {
      apar=ifelse(data_type=="cont", 0.1, 0.1)
      apar=c(apar, 0)
      rnull=function() rnorm(n)
      dnull=function(x) dnorm(x)
      pnull=function(x) pnorm(x)
      palt = function(x, a=apar) {
        if(a == 0)
          pnorm(x)
        else
          pt(x, df=1/a)
      }
      
      ralt = function(a=apar) {
        if(a == 0)
          rnorm(n)
        else
          rt(n, df=1/a)
      }
      TSextra=list(qnull=function(x) qnorm(x))
      bins=seq(-10, 10, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 33
    "normal.outlier1" = (function() {
      apar=ifelse(data_type=="cont", 0.015, 0.085)
      apar=c(apar, 0)
      rnull=function() rnorm(n)
      dnull=function(x) dnorm(x)
      pnull=function(x) pnorm(x)
      palt=function(x, a=apar[1]) (1-a)*pnorm(x)+a*pnorm(x, 3)
      ralt=function(a=apar[1]) 
         c(rnorm(floor((1-a)*n)), rnorm(ceiling(a*n), 3))
      TSextra=list(qnull=function(x) qnorm(x))
      bins=seq(-3, 5, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf 
      collect_case(environment())
    })(),
    # 34
    "normal.outlier2" = (function() {
      apar=ifelse(data_type=="cont", 0.012, 0.03)
      apar=c(apar, 0)
      rnull=function() rnorm(n)
      dnull=function(x) dnorm(x)
      pnull=function(x) pnorm(x)
      palt=function(x, a=apar[1]) 
         (1-2*a)*pnorm(x)+a*punif(x, -5, -3)++a*punif(x, 3, 5)
      ralt=function(a=apar[1]) {
        n1=floor(a*n)
        n2=ceiling(a*n)
        n0=n-n1-n2
        c(rnorm(n0), runif(n1, -5, -3), runif(n2, 3, 5))
      }
      TSextra=list(qnull=function(x) qnorm(x))
      bins=seq(-5, 5, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf   
      collect_case(environment())
    })(),
    # 35
    "normal.normalmixture" = (function() {
      apar=ifelse(data_type=="cont", 0.62, 0.68)
      apar=c(apar, 0)
      rnull=function() rnorm(n)
      dnull=function(x) dnorm(x)
      pnull=function(x) pnorm(x)
      palt=function(x, a=apar[1]) (pnorm(x,-a)+pnorm(x, a))/2
      ralt=function(a=apar[1]) {k=msizes(n, c(1/2, 1/2)); c(rnorm(k[1], -a), rnorm(k[2], a))}
      TSextra=list(qnull=function(x) qnorm(x))
      bins=seq(-3, 3, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf   
      collect_case(environment())
    })(),
    # 36
    "exponential.gamma" = (function() {
      apar=ifelse(data_type=="cont", 1.22, 1.28)
      apar=c(apar, 1)
      rnull=function() rexp(n, 1)
      dnull=function(x) dexp(x, 1)
      pnull=function(x) pexp(x, 1)
      palt=function(x, a=apar[1]) pgamma(x, 1, a)
      ralt=function(a=apar[1]) rgamma(n, 1, a)
      TSextra=list(qnull=function(x) qexp(x, 1))
      bins=seq(0, 10, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 37
    "exponential.weibull" = (function() {
      apar=ifelse(data_type=="cont", 1.22, 1.3)
      apar=c(apar, 1)
      rnull=function() rexp(n, 1)
      dnull=function(x) dexp(x, 1)
      pnull=function(x) pexp(x, 1)
      palt=function(x, a=apar[1]) pweibull(x, 1, a)
      ralt=function(a=apar[1]) rweibull(n, 1, a)
      TSextra=list(qnull=function(x) qexp(x, 1))
      bins=seq(0, 10, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 38
    "exponential.bump" = (function() {
      apar=c(0.124, 0)
      rnull=function() rexp(n, 1)
      dnull=function(x) dexp(x, 1)
      pnull=function(x) pexp(x, 1)
      palt=function(x, a=apar[1]) (1-a)*pexp(x, 1)+a*pnorm(x, 0.5, 0.05)
      ralt=function(a=apar[1]) c(rexp(floor((1-a)*n), 1), rnorm(ceiling(a*n), 0.5, 0.05))
      TSextra=list(qnull=function(x) qexp(x, 1))
      bins=seq(0, 10, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 39
    "chisquare.noncentral" = (function() {
      apar=ifelse(data_type=="cont", 0.6, 3)
      apar=c(apar, 0)
      rnull=function() rchisq(n, 5)
      dnull=function(x) dchisq(x, 5)
      pnull=function(x) pchisq(x, 5)
      palt=function(x, a=apar[1]) pchisq(x, 5, a)
      ralt=function(a=apar[1]) rchisq(n, 5, a)
      TSextra=list(qnull=function(x) qchisq(x, 5))
      bins=seq(0, 25, length=nbins+1)
      bins[nbins+1]=Inf
        
      collect_case(environment())
    })(),
    # 40
    "gamma.gammamixture" = (function() {
      apar=c(0.275, 0)
      rnull=function() rgamma(n, 2, 3)
      dnull=function(x) dgamma(x, 2, 3)
      pnull=function(x) pgamma(x, 2, 3)
      palt=function(x, a=apar[1]) (1-a)*pgamma(x, 2, 3)+a*pgamma(x, 3, 3)
      ralt=function(a=apar[1]) c(rgamma(floor((1-a)*n), 2, 3), rgamma(ceiling(a*n), 3, 3))
      TSextra=list(qnull=function(x) qgamma(x, 2, 3))
      bins=seq(0, 3, length=nbins+1)
      bins[nbins+1]=Inf   
      collect_case(environment())
    })(),
    # 41
    "normal.t.est" = (function() {
      apar=c(0.15, 0)
      rnull = function(p=c(0,1)) rnorm(n, p[1], p[2])
      dnull = function(x, p=c(0,1)) dnorm(x, p[1], p[2])
      pnull = function(x, p=c(0,1)) pnorm(x, p[1], p[2])
      palt = function(x, a=apar) {
        if(a == 0)
          pnorm(x)
        else
          pt(x, df=1/a)
      }
      
      ralt = function(a=apar) {
        if(a == 0)
          rnorm(n)
        else
          rt(n, df=1/a)
      }
      phat = function(x) c(mean(x), sd(x))
      TSextra=list(qnull=function(x, p=c(0,1)) qnorm(x, p[1], p[2]))
      bins=seq(-5, 5, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf
      pnullD=function(p=c(0,1)) pnull(bins[-1],p)
      rnullD=function(p=c(0,1)) 
         c(rmultinom(1, n, diff(pnull(bins, p))))
        
      collect_case(environment())
    })(),
    # 42
    "exponential.weibull.est" = (function() {
      apar=c(1.2, 1)
      pnull = function(x, p=1) pexp(x, p)
      dnull = function(x, p=1) dexp(x, p)
      rnull = function(p=1) rexp(n, p)
      palt=function(x, a=apar[1]) pweibull(x, a, 1)
      ralt = function(a=apar[1]) rweibull(n, a, 1)
      phat = function(x) 1/mean(x)
      TSextra=list(qnull = function(x, p=1) qexp(x, p))
      bins=seq(0, 10, length=nbins+1)
      bins[nbins+1]=Inf
      pnullD=function(p=1) pnull(bins[-1],p)
      rnullD=function(p=1) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      bins=seq(0, 3.5, length=nbins+1)
      bins[nbins+1]=Inf   
      collect_case(environment())
    })(),
    # 43
    "trunc.exponential.linear.est" = (function() {
      apar=c(0.65, 0)
      pnull = function(x, p=1) {
        y=(1-exp(-x*p))/(1-exp(-p))
        ifelse(x<0, 0, ifelse(x>1, 1, y))
        }
      dnull = function(x, p=1) {
        y=p*exp(-x*p)/(1-exp(-p))
        ifelse(x<0, 0, ifelse(x>1, 0, y))
        }
      TSextra=list(qnull = function(x, p=1) -log(1 - x*(1-exp(-p)))/p)
      rnull = function(p=1) {
        x = rexp(2*n/pexp(1, p), p)
        x=x[x<1]
        x[1:n] 
      }
      palt=function(x, a=apar[1]) {
        slope=-1
        p1=ifelse(x<0, 0, ifelse(x>1, 1,(1-slope)*x+slope*x^2))
        y=(1-exp(-x))/(1-exp(-1))
        p2=ifelse(x<0, 0, ifelse(x>1, 1, y))
        (1-a)*p2+a*p1
        }
      ralt = function(a=apar[1]) {
        slope=-1
        y1=(slope-1+sqrt((1-slope)^2+4*slope* runif(n)))/2/slope
        y2 = rexp(2*n/pexp(1, 1), 1)
        y2=y2[y2<1]
        c(y2[1:floor(((1-a)*n))], y1[1:ceiling(a*n)])
      }
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
      }
      bins=seq(0, 1, length=nbins+1)
      pnullD=function(p=1) pnull(bins[-1],p)
      rnullD=function(p=1) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
        
      collect_case(environment())
    })(),
    # 44
    "exponential.gamma.est" = (function() {
      apar=c(1.2, 1)
      pnull = function(x, p=1) pexp(x, p)
      dnull = function(x, p=1) dexp(x, p)
      TSextra=list(qnull = function(x, p=1) qexp(x, p))
      rnull = function(p=1) rexp(n, p)
      palt=function(x, a=apar[1]) pgamma(x, a, 1)
      ralt = function(a=apar[1]) rgamma(n, a, 1)
      phat = function(x) 1/mean(x)
      bins=seq(0, 10, length=nbins+1)
      bins[nbins+1]=Inf
      pnullD=function(p=1) pnull(bins[-1],p)
      rnullD=function(p=1) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
        
      collect_case(environment())
    })(),
    # 45
    "normal.cauchy.est" = (function() {
      apar=c(1, 0)
      xcut=ifelse(data_type=="cont", 3, 3.25)
      pnull = function(x, s=1) {
        y=(pnorm(x, 0, s)-pnorm(-xcut, 0, s))/(2*pnorm(xcut, 0, s)-1)
        y[x<(-xcut)]=0
        y[x>xcut]=1
        y
      }
      dnull = function(x, s=1) {
        y=dnorm(x, 0, s)/(2*pnorm(xcut, 0, s)-1)
        y[x<(-xcut)]=0
        y[x>xcut]=0
        y
      }
      TSextra = list(
        qnull = function(x, s=1) {
          lo <- pnorm(-xcut, 0, s)
          hi <- pnorm( xcut, 0, s)
          qnorm(lo + x*(hi-lo), 0, s)
        }
      )
      rnull = function(s=1) {
        x=rnorm(1.2*n, 0, s)
        x=x[abs(x)<xcut]
        x[1:n]
      }
      # Contamination alternative: a=0 is exactly the truncated N(0,1)
      # member of the composite null; a=1 is truncated standard Cauchy.
      palt=function(x, a=apar[1]) {
         pn=(pnorm(x)-pnorm(-xcut))/(2*pnorm(xcut)-1)
         pc=(pcauchy(x)-pcauchy(-xcut))/(2*pcauchy(xcut)-1)
         y=(1-a)*pn+a*pc
         y[x<(-xcut)]=0
         y[x>xcut]=1
         y
      }
      ralt = function(a=apar[1]) {
         z=numeric(0)
         while(length(z)<n) {
            m=2*(n-length(z))
            x=ifelse(runif(m)<a, rcauchy(m), rnorm(m))
            z=c(z,x[abs(x)<xcut])
         }
         z[1:n]
      }
      bins=seq(-xcut, xcut, length=nbins+1)
      phat = function(x) {
        ll=function(s, x) -sum(log(dnorm(x,0,s))-log(2*pnorm(xcut,0,s)-1))
        optimize(ll, c(0,20), x=x)$minimum
      }
      pnullD=function(s=1) pnull(bins[-1],s)
      rnullD=function(s=1) 
        c(rmultinom(1, n, diff(pnull(bins, s))))
        
      collect_case(environment())
    })(),
    # 46
    "betaa1.betaab.est" = (function() {
      apar=c(1.2, 1)
      pnull = function(x, s=1) pbeta(x, s, 1)
      dnull = function(x, s=1) dbeta(x, s, 1)
      TSextra=list(qnull = function(x, s=1) qbeta(x, s, 1))
      rnull = function(s=1) rbeta(n, s, 1)
      palt=function(x, a=apar[1]) pbeta(x, 1, a)
      ralt = function(a=apar[1]) rbeta(n, 1, a)
      phat = function(x) {
        ol=mean(x)/(1-mean(x))
        t=sum(log(x))
        n=length(x)
        repeat {
          nw=2*ol+t*ol^2/n
          if(abs(ol-nw)<1e-3) break
          ol=nw
        }
        nw
      } 
      bins=seq(0, 1, length=nbins+1)
      pnullD=function(s=1) pnull(bins[-1], s)
      rnullD=function(s=1) 
        c(rmultinom(1, n, diff(pnull(bins, s))))
      
        
      collect_case(environment())
    })(),
    # 47
    "betaa1.betaaa.est" = (function() {
      apar=c(1.2, 1)
      pnull = function(x, s=1) pbeta(x, s, 1)
      dnull = function(x, s=1) dbeta(x, s, 1)
      TSextra=list(qnull = function(x, s=1) qbeta(x, s, 1))
      rnull = function(s=1) rbeta(n, s, 1)
      palt=function(x, a=apar[1]) pbeta(x, a, a)
      ralt = function(a=apar[1]) rbeta(n, a, a)
      phat = function(x) -length(x)/sum(log(x))
      pnullD=function(s=1) pnull(bins[-1], s)
      rnullD=function(s=1) 
        c(rmultinom(1, n, diff(pnull(bins, s))))
      bins=seq(0, 1, length=nbins+1)
      collect_case(environment())
    })(),
    # 48
    "trans-beta.est" = (function() {
      apar=c(2.75, 1)
      #g(x)=(exp(x)-1)/(e-1);g(X)~Beta(a,1)
      dnull=function(x, theta=1) {
        dens <- theta * exp(x) / (exp(1) - 1) *
          ((exp(x) - 1) / (exp(1) - 1))^(theta - 1)
        dens[x < 0 | x > 1] <- 0
        dens
      }
      pnull=function(x, theta=1) {
        p <- ((exp(x) - 1) / (exp(1) - 1))^theta
        p[x <= 0] <- 0
        p[x >= 1] <- 1
        p
      }
      
      rnull=function(theta=1) {
        log(1+(exp(1)-1)*runif(n)^(1/theta))
      }  
      palt <- function(x, a=apar[1]) {
        ifelse(x <= 0, 0,
               ifelse(x >= 1, 1,
                      ((exp(a * x) - 1) / (exp(a) - 1))))
      }
      ralt <- function(a=apar[1]) 
        log(1 + (exp(a)-1)*runif(n))/a
      
      TSextra=list(qnull=function(x, theta=1) {
        log(1+(exp(1) - 1)*x^(1/theta))
      })
      phat <- function(x) {
        -length(x)/sum(log((exp(x)-1)/(exp(1)-1)))
      }
      pnullD=function(p=1) pnull(bins[-1], p)
      rnullD=function(p=1) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      bins=seq(0, 1, length=nbins+1)

      collect_case(environment())
    })(),
    # 49
    "double.exponential.est" = (function() {
      apar=ifelse(data_type=="cont", 1.3, 1.4)
      apar=c(apar, 1)
      rnull=function(a=1) ifelse(runif(n)<0.5,-rexp(n, a), rexp(n, a))
      dnull=function(x, a=1) a*exp(-a*abs(x))/2
      pnull=function(x, a=1) 
         ifelse(x<0, exp(a*x)/2, 1-exp(-a*x)/2)
      palt=function(x, a=apar[1]) 
        ifelse(x<0, (1-pgamma(-x, a, 1))/2, (1+pgamma(x, a, 1))/2)
      ralt=function(a=apar[1]) 
        ifelse(runif(n)<1/2, -rgamma(n, a, 1), rgamma(n, a, 1))
      TSextra=list(qnull=function(x, a=1) {
        x[x<0.00001]=0.00001
        x[x>1-0.00001]=1-0.00001
        ifelse(x<1/2, log(2*x)/a, -log(2*(1-x))/a)
      })
      phat=function(x) 1/mean(abs(x))
      bins=seq(-10, 10, length=nbins+1)
      bins[1]=-Inf;bins[nbins+1]=Inf
      pnullD=function(p=1) pnull(bins[-1],p)
      rnullD=function(p=1) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      
        
      collect_case(environment())
    })(),
    # 50
    "pareto.est" = (function() {
      apar=c(0.8, 0)
      dnull <- function(x, alpha=10) {
        ifelse(x<0, 0, alpha/(x+1)^(alpha + 1))
        }
      pnull <- function(x, alpha=10) ifelse(x<0, 0, 1-(1/(x+1))^alpha)
      TSextra=list(qnull = function(x, alpha=10) {
        (1 - x)^(-1 / alpha)-1})
      rnull <- function(alpha=10) runif(n)^(-1 / alpha)-1
      # Contamination alternative: a=0 is exactly Pareto(10); a=1 is
      # the original Gamma(shape=1.07, rate=8) departure.
      palt=function(x, a=apar[1]) (1-a)*pnull(x,10)+a*pgamma(x,1.07,8)
      ralt=function(a=apar[1]) {
         I=runif(n)<a
         x=numeric(n)
         x[!I]=runif(sum(!I))^(-1/10)-1
         x[I]=rgamma(sum(I),1.07,8)
         x
      }
      phat=function(x) length(x)/sum(log(x+1))
      bins=seq(0, 1, length=nbins+1)
      bins[nbins+1]=Inf
      pnullD=function(p=10) pnull(bins[-1],p)
      rnullD=function(p=10) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      
        
      collect_case(environment())
    })(),
    # 51
    "betamixture1.est" = (function() {
      apar=c(2.5, 2)
      pnull = function(x, a=1/2) a*pbeta(x, 1, 1)+(1-a)*pbeta(x, 2, 2)
      dnull = function(x, a=1/2) a*dbeta(x, 1, 1)+(1-a)*dbeta(x, 2, 2)
      TSextra=list(qnull = function(x, a=1/2) finv(x, function(y) pnull(y, a), c(0, 1)))
      rnull = function(a=1/2) c(rbeta(floor(a*n), 1, 1), rbeta(ceiling((1-a)*n), 2, 2))
      palt=function(x, a=3) pbeta(x, 1, 1)/2+pbeta(x, a, a)/2
      ralt = function(a=apar[1]) 
         c(rbeta(floor(n/2), 1, 1), rbeta(ceiling(n/2), a, a))      
      phat = function(x) Rgof::mlemix(dbeta(x, 1 ,1), dbeta(x, 2, 2))
      bins=seq(0, 1, length=nbins+1)
      pnullD=function(p=1/2) pnull(bins[-1],p)
      rnullD=function(p=1/2) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      collect_case(environment())
    })(),
    # 52
    "betamixture2.est" = (function() {
      apar=c(1.15, 0)
      pnull = function(x, a=1/2) a*pbeta(x, 3, 1)+(1-a)*pbeta(x, 1, 3)
      dnull = function(x, a=1/2) a*dbeta(x, 3, 1)+(1-a)*dbeta(x, 1, 3)
      TSextra=list(qnull = function(x, a=1/2) finv(x, function(y) pnull(y, a), c(0, 1)))
      rnull = function(a=1/2) c(rbeta(floor(a*n), 3, 1), rbeta(ceiling((1-a)*n), 1, 3))
      palt=function(x, a=apar[1]) 
        pbeta(x, 3+a, 1)/2+pbeta(x, 1, 3)/2
      ralt = function(a=apar[1]) 
        c(rbeta(floor(n/2), 3+a, 1), rbeta(ceiling(n/2), 1, 3))      
      phat = function(x) Rgof::mlemix(dbeta(x, 3 ,1), dbeta(x, 1, 3))
      bins=seq(0, 1, length=nbins+1)
      pnullD=function(p=1/2) pnull(bins[-1],p)
      rnullD=function(p=1/2) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      collect_case(environment())
    })(),
    # 53
    "betamixture3.est" = (function() {
      apar=c(0.13, 0)
      pnull = function(x, a=1/2) a*pbeta(x, 5, 2)+(1-a)*pbeta(x, 2, 5)
      dnull = function(x, a=1/2) a*dbeta(x, 5, 2)+(1-a)*dbeta(x, 2, 5)
      TSextra=list(qnull = function(x, a=1/2) 
        finv(x, function(y) pnull(y, a), c(0, 1)))
      rnull = function(a=1/2) 
        c(rbeta(floor(a*n), 5, 2), rbeta(ceiling((1-a)*n), 2, 5))
      palt=function(x, a=apar[1]) 
        pbeta(x, 5+a, 2-a)/2+pbeta(x, 2-a, 5+a)/2
      ralt = function(a=apar[1]) 
        c(rbeta(floor(n/2), 5+a, 2-a), rbeta(ceiling(n/2), 2-a, 5+a))      
      phat = function(x) Rgof::mlemix(dbeta(x, 5 ,2), dbeta(x, 2, 5))
      bins=seq(0, 1, length=nbins+1)
      pnullD=function(p=1/2) pnull(bins[-1],p)
      rnullD=function(p=1/2) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      collect_case(environment())
    })(),
    # 54
    "normalmixture.est" = (function() {
      apar=c(2.3, 2)
      pnull = function(x, a=1/2) a*pnorm(x)+(1-a)*pnorm(x, 2)
      dnull = function(x, a=1/2) a*dnorm(x)+(1-a)*dnorm(x, 2)
      TSextra=list(qnull = function(x, a=1/2)  finv(x, function(y) pnull(y, a), c(-10, 10)))
      rnull = function(a=1/2) c(rnorm(floor(a*n)), rnorm(ceiling((1-a)*n), 2))
      palt=function(x, a=apar[1]) 0.5*pnorm(x)+0.5*pnorm(x, a)
      ralt = function(a=apar[1]) {k=msizes(n, c(1/2, 1/2)); c(rnorm(k[1]), rnorm(k[2], a))}
      phat = function(x) Rgof::mlemix(dnorm(x), dnorm(x, 2))
      bins=seq(-3, 5, length=nbins+1)
      bins[1]=-Inf
      bins[nbins+1]=Inf
      pnullD=function(p=1/2) pnull(bins[-1],p)
      rnullD=function(p=1/2) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      collect_case(environment())
    })(),
    #55
    "laplace1.est" = (function() {
      apar=c(0.8, 0)
      dnull=function(x, p=c(0, 1)) 
        exp(-abs(x - p[1]) / p[2]) / (2 * p[2])
      pnull=function(x, p=c(0,1)) {
        ifelse(x < p[1], 0.5 * exp((x - p[1]) / p[2]),
               1 - 0.5 * exp(-(x - p[1]) / p[2]))
      }
      TSextra=list(qnull=function(x, p=c(0,1)) {
        ifelse(x < 0.5, p[1] + p[2] * log(2 * x),
               p[1] - p[2] * log(2 * (1 - x)))
      })
      rnull=function(p=c(0,1))  {
        x=runif(n)
        ifelse(x < 0.5, p[1] + p[2] * log(2 * x),
               p[1] - p[2] * log(2 * (1 - x)))
      }
      phat=function(x) c(stats::median(x), mean(abs(x - stats::median(x))))
      # Contamination alternative: a=0 is exactly Laplace(0,1); a=1 is
      # the fixed three-component normal-mixture departure.
      palt <- function(x, a=apar[1]) {
        sd <- sqrt(qexp(c(.1,.3,.5), rate=1))
        G <- 0.6*pnorm(x,0,sd[1]) +
             0.2*pnorm(x,0,sd[2]) +
             0.2*pnorm(x,0,sd[3])
        (1-a)*pnull(x,c(0,1))+a*G
      }
      
      ralt <- function(a=apar[1]) {
        I <- runif(n)<a
        x <- numeric(n)
        u <- runif(sum(!I))
        x[!I] <- ifelse(u<.5,log(2*u),-log(2*(1-u)))
        m <- sum(I)
        if(m) {
          sd <- sqrt(qexp(c(.1,.3,.5), rate=1))
          k <- sample(1:3,m,replace=TRUE,prob=c(.6,.2,.2))
          x[I] <- rnorm(m,0,sd[k])
        }
        x
      }   
      bins=seq(-5, 5, length=nbins+1)
      bins[1]=-Inf
      bins[nbins+1]=Inf
      pnullD=function(p=c(0,1)) pnull(bins[-1],p)
      rnullD=function(p=c(0, 1)) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      phatD=function(x) grouped_mle_laplace(x, bins, pnull, phat)
      collect_case(environment())
    })(),
    #56
    "laplace2.est" = (function() {
      apar=c(1.2, 1)
      dnull=function(x, p=c(0, 1)) exp(-abs(x - p[1]) / p[2]) / (2 * p[2])
      pnull=function(x, p=c(0,1)) {
        ifelse(x < p[1], 0.5 * exp((x - p[1]) / p[2]),
               1 - 0.5 * exp(-(x - p[1]) / p[2]))
      }
      TSextra=list(qnull=function(x, p=c(0,1)) {
        ifelse(x < 0.5, p[1] + p[2] * log(2 * x),
               p[1] - p[2] * log(2 * (1 - x)))
      })
      rnull=function(p=c(0,1))  {
        x=runif(n)
        ifelse(x < 0.5, p[1] + p[2] * log(2 * x),
               p[1] - p[2] * log(2 * (1 - x)))
      }
      phat=function(x) c(stats::median(x), mean(abs(x - stats::median(x))))
      palt <- function(x, a=apar[1]) {
        ifelse(x<0, 1-pgamma(-x, a, 1),
               1+pgamma(x, a, 1))/2
      }
      ralt=function(a=apar[1]) ifelse(runif(n)<0.5,-1,1)*rgamma(n, a, 1)     
      bins=seq(-5, 5, length=nbins+1)
      bins[1]=-Inf
      bins[nbins+1]=Inf
      pnullD=function(p=c(0, 1)) pnull(bins[-1],p)
      rnullD=function(p=c(0, 1)) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      phatD=function(x) grouped_mle_laplace(x, bins, pnull, phat)
      collect_case(environment())
    })(),
    #57
    "truncexp1.est" = (function() {
      apar=c(1.3, 1)
      pnull=function(x, p=c(1,2)) {
        ifelse(x<p[2], pexp(x, p[1])/pexp(p[2], p[1]), 1) 
      }
      dnull=function(x, p=c(1,2)) 
        ifelse(x<p[2], dexp(x, p[1])/pexp(p[2], p[1]), 0)
      rnull= function(p=c(1,2)) {
        x=rexp(2*n/pexp(p[2], p[1]) , p[1])
        x[x<p[2]][1:n]
      }
      phat=Rgof::mletexp
      TSextra=list(qnull=function(x, p=c(1,2)) 
        -1/p[1]*log(1-x*(1-exp(-p[1]*p[2]))))
      palt <- function(x,a=apar[1])
        ifelse(x < 0, 0,
               ifelse(x > 2, 1,
                      pgamma(x,a,1)/pgamma(2,a,1)))
      ralt=function(a=apar[1]) {
        x=rgamma(2*n/pgamma(2, a, 1), a, 1)
        x[x<2][1:n]
      }
      bins=seq(0, 2, length=nbins+1)
      pnullD=function(p=c(1,1)) pnull(bins[-1],p)
      rnullD=function(p=c(1,1)) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      # grouped-data mle: the truncation point is the upper edge of the last
      # nonempty bin, the rate maximizes the multinomial likelihood given that
      phatD=function(x) {
        up=bins[max(which(x>0))+1]
        nll=function(r) {
          pr=diff(pnull(bins, c(r, up)))
          -sum(x[x>0]*log(pmax(pr[x>0], 1e-300)))
        }
        c(stats::optimize(nll, c(0.01, 20))$minimum, up)
      }
      collect_case(environment())
    })(),
    #58
    "truncexp2.est" = (function() {
      apar=c(0.45, 0)
      pnull=function(x, p=1) {
          ifelse(x<1, 0,
            ifelse(x>2, 1,
                (pexp(x,p)-pexp(1,p))/(pexp(2,p)-pexp(1,p))))
      }  
      dnull=function(x, p=1) 
        ifelse(x<1|x>2, 0 ,dexp(x, p)/(pexp(2, p)-pexp(1, p)))
      rnull= function(p=1) {
        x=rexp(2*n/(pexp(2, p)-pexp(1, p)), p)
        x[x>1 & x<2][1:n]
      }
      phat=Rgof::mledexp
      TSextra=list(qnull=function(x, p=1) 
        -1/p*log( exp(-p) -(exp(-p)-exp(-2*p))*x))
      palt=function(x, a=apar[1]) {
        p1=ifelse(x<1, 0,
           ifelse(x>2, 1,
            (pexp(x, 1)-0.6321206)/0.2325442))
        p2=pbeta(x-1, 2, 1)
        (1-a)*p1+a*p2
      }
      ralt= function(a=apar[1]) {
        x=rexp(floor(8.600517*n), 1)
        x=x[x>1 & x<2]
        c(x[1:ceiling((1-a)*n)], rbeta(floor(a*n), 2, 1)+1)
      }
      bins=seq(1, 2, length=nbins+1)
      pnullD=function(p=1) pnull(bins[-1],p)
      rnullD=function(p=1) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      collect_case(environment())
    })(),
    # 59
    "gammamixture1.est" = (function() {
      apar=c(1.75, 1)
      pnull = function(x, a=2/5) a*pgamma(x, 1, 1)+(1-a)*pgamma(x, 3, 1)
      dnull = function(x, a=2/5) a*dgamma(x, 1, 1)+(1-a)*dgamma(x, 3, 1)
      TSextra=list(qnull = function(x, a=2/5) finv(x, function(y) pnull(y, a), c(0, 50)))
      rnull = function(a=2/5) c(rgamma(floor(a*n), 1 ,1), rgamma(ceiling((1-a)*n), 3, 1))
      palt=function(x, a=apar[1]) 0.4*pgamma(x, 1, a)+0.6*pgamma(x, 3, 1)
      ralt = function(a=apar[1]) c(rgamma(floor(2*n/5), 1 ,a), rgamma(ceiling(3*n/5), 3, 1))      
      phat = function(x) Rgof::mlemix(dgamma(x,1,1), dgamma(x,3,1))
      bins=seq(0, 10, length=nbins+1)
      bins[nbins+1]=Inf
      pnullD=function(p=1/2) pnull(bins[-1],p)
      rnullD=function(p=1/2) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      collect_case(environment())
    })(),
    # 60
    "gammamixture2.est" = (function() {
      apar=c(9.35, 10)
      pnull = function(x, a=0.3) a*pgamma(x, 3, 1)+(1-a)*pgamma(x, 10, 1)
      dnull = function(x, a=0.3) a*dgamma(x, 3, 1)+(1-a)*dgamma(x, 10, 1)
      TSextra=list(qnull = function(x, a=0.3) finv(x, function(y) pnull(y, a), c(0, 50)))
      rnull = function(a=0.3) c(rgamma(floor(a*n), 3 ,1), rgamma(ceiling((1-a)*n), 10, 1))
      palt=function(x, a=apar[1]) 0.3*pgamma(x, 3, 1)+0.7*pgamma(x, a, 1)
      ralt = function(a=apar[1]) c(rgamma(floor(3*n/10), 3 ,1), rgamma(ceiling(7*n/10), a, 1))      
      phat = function(x) Rgof::mlemix(dgamma(x,3,1), dgamma(x,10,1))
      bins=seq(0, 25, length=nbins+1)
      bins[nbins+1]=Inf
      pnullD=function(p=1/2) pnull(bins[-1],p)
      rnullD=function(p=1/2) 
        c(rmultinom(1, n, diff(pnull(bins, p))))
      collect_case(environment())
    })()
  )

  if (is.null(study))
    stop("unknown case study: ", which, call.=FALSE)

  pnull <- study$pnull
  dnull <- study$dnull
  rnull <- study$rnull
  TSextra <- study$TSextra
  if (!is.null(study$palt)) palt <- study$palt
  if (!is.null(study$ralt)) ralt <- study$ralt
  if (!is.null(study$phat)) phat <- study$phat
  phatD <- study$phatD
  if (!is.null(study$pnullD)) pnullD <- study$pnullD
  if (!is.null(study$rnullD)) rnullD <- study$rnullD
  if (!is.null(study$bins)) bins <- study$bins
  if (!is.null(study$Ra)) Ra <- study$Ra
  if (!is.null(study$apar)) apar <- study$apar

# End of definitions  
#   

# Create output 

  # Keep the alternative parameter explicit, as in R2sample::case.studies().
  # ralt() and palt() both use apar as their default alternative parameter.
  ralt1=function() ralt(apar[1])
  palt1=function(x) palt(x, apar[1])

  Ra=bins[c(1, nbins+1)]
  if(data_type=="disc") {
    vals=(bins[-1]+bins[-(nbins+1)])/2 # midpoints of intervals
    phat1=function(x) {
      y=rep(vals, x)
      phat(y)
    }
    # Where the continuous estimator applied to the bin midpoints is lattice valued
    # (median, maximum), and the parametric bootstrap therefore fails, the case
    # study has the grouped-data mle phatD for discrete data
    if(!is.null(phatD)) phat1=phatD
    if(endsWith(which,"est")) {
      pnull1=pnullD
      rnull=rnullD
    }
    else {
      pnull1=function() pnull(bins[-1])
      rnull=function() c(rmultinom(1, n, diff(pnull(bins))))
    }
    vals = (bins[-1]+bins[-(nbins+1)])/2
    if(is.infinite(vals[1])) vals[1]=bins[2]-1
    if(is.infinite(vals[length(vals)]))
      vals[length(vals)]=bins[length(bins)-1]+1

    # For discrete data ralt returns cell counts, but still accepts the
    # alternative parameter with default apar.
    raltD=function(a=apar[1]) c(rmultinom(1, n, diff(palt(bins, a))))

    return(list(pnull=pnull1, rnull=rnull,
                palt=palt, ralt=raltD, apar=apar[1], apar0=apar[2],
                phat=phat1, ralt1=raltD, vals=vals,
                TSextra=TSextra, Range=Ra))
  }
  return(list(pnull=pnull, rnull=rnull, TSextra=TSextra,
              phat=phat, dnull=dnull, Range=Ra,
              palt=palt, ralt=ralt, apar=apar[1], apar0=apar[2],
              ralt1=ralt1, vals=NA))
}

# grouped-data mle of the Laplace distribution (location, scale), for binned data:
# maximizes the multinomial log-likelihood sum x_j log(F(b_j)-F(b_{j-1})), starting
# at the estimates of phat applied to the bin midpoints.
grouped_mle_laplace <- function(x, bins, pnull, phat) {
  k <- length(bins)-1
  mids <- (bins[-1]+bins[-(k+1)])/2
  if(is.infinite(mids[1])) mids[1] <- bins[2]-1
  if(is.infinite(mids[k])) mids[k] <- bins[k]+1
  st <- phat(rep(mids, x))
  nll <- function(q) {
    pr <- diff(pnull(bins, c(q[1], exp(q[2]))))
    -sum(x[x>0]*log(pmax(pr[x>0], 1e-300)))
  }
  o <- stats::optim(c(st[1], log(max(st[2], 0.01))), nll)
  c(o$par[1], exp(o$par[2]))
}

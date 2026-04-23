####packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","data.table","latex2exp", "Cairo", "gamlss")


#-------------------------------------------------------------------------------
#expression of log density function 
l_DTED <- expression(log(
((sigma+log(2))/(mu))*exp(-(y*(sigma+log(2)))/(mu))*(1- (1/(sigma+log(2)))*log(log((1-exp(-(sigma+log(2))))/0.5))*(1-exp(-((sigma+log(2))*y)/mu))*((log((1-exp(-(sigma+log(2))))/(0.5)))^(y/mu))*exp((y*(sigma+log(2)))/(mu)))*exp(-((log((1-exp(-(sigma+log(2))))/0.5))^(y/mu)))
))

######first-order derivate with respect to each parameter
m1 <- D(l_DTED, "mu")
s1 <- D(l_DTED, "sigma")

#second-order derivate
#ms2 <- D(m1, "sigma") #igual a sm2

m2  <- D(m1, "mu")
s2  <- D(s1, "sigma")
ms2 <- D(m1, "sigma")


#------------------------------------------------------------------------------
# Dutta-exponential distribution in GAMLSS
#

DTED <- function(mu.link = "log", sigma.link = "log"){
  
  mstats <- checklink("mu.link", "DTED", substitute(mu.link),
                      c("log","inverse","identity", "sqrt", "own"))
  
  dstats <- checklink("sigma.link", "DTED", substitute(sigma.link),
                      c("log", "inverse", "identity", "sqrt", "logit", "own"))
    
  structure(
    list(family = c("DTED", " Dutta-exponential distribution"),
         parameters = list(mu=TRUE, sigma=TRUE),
         nopar = 2, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
                  
         ######derivates
         dldm =  function(y, mu, sigma){
           dldm = eval(m1)
           dldm
         },
         dldd = function(y, mu, sigma){
           dldd = eval(s1)
           dldd
         },
                

        d2ldm2 = function(y, mu, sigma) {
          d2ldm2 <- eval(m2)
          d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
          d2ldm2
        },

        d2ldd2 = function(y, mu, sigma) {
          d2ldd2 <- eval(s2)
          d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
          d2ldd2
        },

        d2ldmdd = function(y, mu, sigma) {
          d2ldmdd <- eval(ms2)
          d2ldmdd[is.na(d2ldmdd)] <- 0
          d2ldmdd
        },
  
           
         #####
         G.dev.incr = function(y,mu,sigma, ...) -2*dDTED(y = y, mu = mu, sigma = sigma , log=TRUE),
         
         rqres = expression(rqres(pfun="pDTED", type = "Continuous", y=y, mu=mu, sigma=sigma)),
         
         #####Initial values for mu and sigma
         mu.initial = expression(mu <-  rep(median(y),length(y))),
         sigma.initial = expression(sigma <- rep(1.0,length(y))),
                
         #####Restrictions
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

# Density function
dDTED <- function(y, mu = 1, sigma = 1, log = FALSE) {
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(sigma <= 0)) stop("sigma must be positive")
  if (any(y <= 0)) stop("y must be positive")

  A <- sigma + log(2)
  B <- log((1 - exp(-A)) / 0.5)
  C <- log(B)

  fy1 <- (A/mu) *
    exp(-(y*A)/mu) *
    (
      1 -
      (1/A) *
      C *
      (1 - exp(-(A*y)/mu)) *
      (B^(y/mu)) *
      exp((y*A)/mu)
    ) *
    exp(-(B^(y/mu)))

  if (log) {
    fy <- log(fy1)
    fy <- ifelse(y <= 0, -Inf, fy)
  } else {
    fy <- ifelse(y <= 0, 0, fy1)
  }

  return(fy)
}

#dDTED(y=1, mu=1, sigma=1)

#-------------------------------------------------------------------------------
# Cumulative function
pDTED <- function(q, mu=1, sigma=1, lower.tail = TRUE, log = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q <= 0)) stop(paste("y must be positive", "\n", ""))
  
  cdf <- (1-exp(-(q*(sigma+log(2)))/mu))*exp(-((log((1-exp(-(sigma+log(2))))/.5))^(q/mu)))

  if(lower.tail==TRUE) cdf <- cdf else cdf <- 1 - cdf
  if(log==TRUE) cdf <- log(cdf)
  
  return(cdf)
}

#pDTED(0.5, mu=1, sigma=1)

#-------------------------------------------------------------------------------
# Quantile function
qDTED <- function(p, mu=1, sigma=1, lower.tail = TRUE, log = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  
  if (log==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  root=NULL
  for(i in 1:length(p)){
    
    prob=p[i]
    
    f=function(y, mu, sigma, p) pDTED(y, mu, sigma) - p
    
    root[i]=uniroot(f,c(0.001, 100000), tol = 0.0001, mu, sigma, p)$root
  }
  return(root)
}


# qDTED(0.5, mu=5, sigma=1)


#-------------------------------------------------------------------------------
# Pseudo-random function
rDTED <- function(n, mu = 1, sigma = 1){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  
  y <- numeric(0)
  
  for(i in 1:n) 
  {
    u <-  runif(1)
    f <- function(y, sigma, mu, u) pDTED(y, mu, sigma) - u
    
    y[i] <- uniroot(f,c(0.001, 100000), tol = 0.0001, sigma, mu, u)$root
  }
  
  return(y)
}


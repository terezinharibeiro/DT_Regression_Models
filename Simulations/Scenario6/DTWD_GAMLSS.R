####packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","data.table","latex2exp", "Cairo", "gamlss")

#-------------------------------------------------------------------------------
#expression of log density function 
l_DTWD <- expression(log(
nu * (((sigma + (log(2))^(1/nu)) / mu)^nu) * (y^(nu - 1)) *
exp(-(((sigma + (log(2))^(1/nu)) * y / mu)^nu)) * (1 - (log(log((1 - exp(-((sigma + (log(2))^(1/nu))^nu))) / 0.5)) / ((sigma + (log(2))^(1/nu))^nu)) *
(1 - exp(-(((sigma + (log(2))^(1/nu)) * y / mu)^nu))) * ((log((1 - exp(-((sigma + (log(2))^(1/nu))^nu))) / 0.5))^((y / mu)^nu)) *  exp(((sigma + (log(2))^(1/nu)) * y / mu)^nu)) * exp(-(log((1 - exp(-((sigma + (log(2))^(1/nu))^nu))) / 0.5))^((y / mu)^nu))
))

#------------------------------------------------------------------------------
# Dutta-Weibull distribution in GAMLSS
#

DTWD <- function(mu.link = "log", sigma.link = "log", nu.link = "log"){
  
  mstats <- checklink("mu.link", "DTWD", substitute(mu.link),
                      c("log","inverse","identity", "sqrt", "own"))
  
  dstats <- checklink("sigma.link", "DTWD", substitute(sigma.link),
                      c("log", "inverse", "identity", "sqrt", "logit", "own"))

  vstats <- checklink("nu.link", "DTWD", substitute(nu.link),
                      c("log","inverse","identity", "sqrt", "own"))
    
  structure(
    list(family = c("DTWD", " Dutta-Weibull distribution"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE),
         nopar = 3, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         nu.link = as.character(substitute(nu.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         nu.dr = vstats$mu.eta,

         
         ######derivates
         dldm =  function(y, mu, sigma, nu){
           m1 <- D(l_DTWD, "mu")
           dldm = eval(m1)
           dldm
         },
         dldd = function(y, mu, sigma, nu){
           s1 <- D(l_DTWD, "sigma")
           dldd = eval(s1)
           dldd
         },
         dldv = function(y, mu, sigma, nu){
           n1 <- D(l_DTWD, "nu")
           dldv = eval(n1)
           dldv
         },       

        d2ldm2 = function(y, mu, sigma, nu) {
          m1 <- D(l_DTWD, "mu")
          m2  <- D(m1, "mu")
          d2ldm2 <- eval(m2)
          d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
          d2ldm2
        },

        d2ldd2 = function(y, mu, sigma, nu) {
          s1 <- D(l_DTWD, "sigma")
          s2  <- D(s1, "sigma")
          d2ldd2 <- eval(s2)
          d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
          d2ldd2
        },
        
        d2ldv2 = function(y, mu, sigma, nu) {
          n1 <- D(l_DTWD, "nu")
          n2  <- D(n1, "nu")
          d2ldv2 <- eval(n2)
          d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2, -1e-15)
          d2ldv2
        },
        
        d2ldmdd = function(y, mu, sigma, nu) {
          m1 <- D(l_DTWD, "mu")
          ms2 <- D(m1, "sigma")
          d2ldmdd <- eval(ms2)
          d2ldmdd[is.na(d2ldmdd)] <- 0
          d2ldmdd
        },
  
       d2ldmdv = function(y, mu, sigma, nu) {
          m1 <- D(l_DTWD, "mu")
          mn2 <- D(m1, "nu")
          d2ldmdv <- eval(mn2)
          d2ldmdv[is.na(d2ldmdv)] <- 0
          d2ldmdv
        },

       d2ldddv =  function(y, mu, sigma, nu) {
          s1 <- D(l_DTWD, "sigma")
          sn2 <- D(s1, "nu")
          d2ldddv <- eval(sn2)
          d2ldddv[is.na(d2ldddv)] <- 0
          d2ldddv
        },



    
         #####
         G.dev.incr = function(y,mu,sigma,nu, ...) -2*dDTWD(y = y, mu = mu, sigma = sigma, nu=nu, log=TRUE),
         
         rqres = expression(rqres(pfun="pDTWD", type = "Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
         
         #####Initial values for mu, sigma and nu
         mu.initial = expression(mu <-  rep(median(y),length(y))),
         sigma.initial = expression(sigma <- rep(1.0,length(y))),
         nu.initial = expression(nu <- rep(1.0,length(y))),
       
         #####Restrictions
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu) all(nu > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

# Density function
dDTWD <- function(y, mu = 1, sigma = 1, nu=1, log = FALSE) {
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(sigma <= 0)) stop("sigma must be positive")
  if (any(nu <= 0)) stop("nu must be positive")
  if (any(y <= 0)) stop("y must be positive")

  c <- sigma + (log(2))^(1/nu)
  A <- log((1 - exp(-(c^nu))) / 0.5)
  t1 <- (c * y / mu)^nu
  t2 <- (y / mu)^nu
  
  fy1 <- nu * ((c / mu)^nu) * (y^(nu - 1)) *
    exp(-t1) * (1 - (log(A) / (c^nu)) *
       (1 - exp(-t1)) * (A^t2) *  exp(t1)) * exp(-A^t2)

  if (log) {
    fy <- log(fy1)
    fy <- ifelse(y <= 0, -Inf, fy)
  } else {
    fy <- ifelse(y <= 0, 0, fy1)
  }

  return(fy)
}

#dDTWD(y=1, mu=1, sigma=1, nu=1)

#-------------------------------------------------------------------------------
# Cumulative function
pDTWD <- function(q, mu=1, sigma=1, nu=1, lower.tail = TRUE, log = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) stop(paste("nu must be positive", "\n", ""))
  if (any(q <= 0)) stop(paste("y must be positive", "\n", ""))
  
  c <- sigma + (log(2))^(1/nu)
  A <- log((1 - exp(-(c^nu))) / 0.5)
  t1 <- (c * q / mu)^nu
  t2 <- (q / mu)^nu

  cdf <- exp(-A^t2)*(1-exp(-t1))

  if(lower.tail==TRUE) cdf <- cdf else cdf <- 1 - cdf
  if(log==TRUE) cdf <- log(cdf)
  
  return(cdf)
}

#pDTWD(1, mu=1, sigma=1, nu=1)

#-------------------------------------------------------------------------------
# Quantile function
qDTWD <- function(p, mu=1, sigma=1, nu=1, lower.tail = TRUE, log = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) stop(paste("nu must be positive", "\n", ""))

  if (log==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  root=NULL
  for(i in 1:length(p)){
    
    prob=p[i]
    
    f=function(y, mu, sigma, nu, p) pDTWD(y, mu, sigma, nu) - p
    
    root[i]=uniroot(f,c(0.001, 100000), tol = 0.0001, mu, sigma, nu, p)$root
  }
  return(root)
}


# qDTWD(0.5, mu=5, sigma=1, nu=1)


#-------------------------------------------------------------------------------
# Pseudo-random function
rDTWD <- function(n, mu = 1, sigma = 1, nu=1){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) stop(paste("nu must be positive", "\n", ""))
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  
  y <- numeric(0)
  
  for(i in 1:n) 
  {
    u <-  runif(1)
    f <- function(y, mu, sigma, nu, u) pDTWD(y, mu, sigma, nu) - u
    
    y[i] <- uniroot(f,c(0.001, 100000), tol = 0.0001, mu, sigma, nu, u)$root
  }
  
  return(y)
}

# rDTWD(10, mu=5, sigma=1, nu=1)

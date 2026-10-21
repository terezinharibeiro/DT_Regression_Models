
#**************************************************************************************************************************************************************************************#

FunSample2 <- function(n, X, Z, W, theta, linkmu, linksigma, linknu){
   
  kk1 <- ncol(X); kk2 <- ncol(Z); kk3 <- ncol(W)
  
  # vector of coefficient regressors
  alpha <- theta[1:kk1]
  gama <- theta[(kk1+1.0):(kk1+kk2)]
  xi <- theta[(kk1+kk2+1.0):(kk1+kk2+kk3)]
  
  #linear predictors                                                
  eta <- as.vector(X%*%alpha)	
  delta <- as.vector(Z%*%gama)
  zeta <-  as.vector(W%*%xi) 
  
  if(linkmu == "log") mu <- exp(eta) #log
  if(linkmu == "identity") mu <- eta #identity
  if(linkmu == "sqrt") mu <- eta^2 #sqrt
  if(linkmu == "inverse") mu <- 1/eta #inverse
  
  if(linksigma == "log") sigma <- exp(delta) #log
  if(linksigma == "identity") sigma <- delta #identity
  if(linksigma == "sqrt") sigma <- delta^2 #sqrt
  if(linksigma == "inverse") sigma <-  1/delta #inverse
  
  
  if(linknu == "log") nu <- exp(zeta) #log
  if(linknu == "identity") nu <- zeta #identity
  if(linknu == "sqrt") nu <- zeta^2 #sqrt
  if(linknu == "inverse") nu <-  1/zeta #inverse
  

 
 #Generating y
  y <- numeric(0)
  
  for(i in 1:n){
    y[i] <- rDTWD(1, mu[i], sigma[i], nu[i])                                            
  }
  
  results <- list(y=y, X=X, Z=Z, W=W, mu=mu, sigma= sigma, nu=nu)
  return(results)
}
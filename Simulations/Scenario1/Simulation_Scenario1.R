####pacotes
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load("tidyverse","data.table","latex2exp", "Cairo", "gamlss", 
 #              "devtools", "remotes")

search()


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

source("DTWD_GAMLSS.r")
source("SAMPLE2.r")


####txt file names
#
sink(str_c("Results_Scenario", str_sub(getwd(), -1, -1), ".txt"))

#
names_results1 <- c(str_c("Estimates_Scenario_", str_sub(getwd(), -1, -1), "_n50"), 
                    str_c("Estimates_Scenario_", str_sub(getwd(), -1, -1), "_n100"), 
                    str_c("Estimates_Scenario_", str_sub(getwd(), -1, -1), "_n300"))

#
names_results2 <- c(str_c("SE_Estimates_Scenario_", str_sub(getwd(), -1, -1), "_n50"), 
                    str_c("SE_Estimates_Scenario_", str_sub(getwd(), -1, -1), "_n100"), 
                    str_c("SE_Estimates_Scenario_", str_sub(getwd(), -1, -1), "_n300"))

#-------------------------------Global Variables-------------------------------#
VN <- c(50, 100, 300) #Sample sizes
VALPHA <- c(1.4, 1.0) #true alpha values - MU
VGAMA <- c(-2.5) #true gamma values - SIGMA
VXI <- c(1.0) #true xi values - NU
VTHETA <- c(VALPHA, VGAMA, VXI) # true theta values
NREP <- 1000 #number of Monte Carlo replicates 

kk1 <- length(VALPHA)
kk2 <- length(VGAMA) 
kk3 <- length(VXI) 


#--------------------------------begin loop------------------------------------#

for(l in 1:length(VN)){
  N = VN[l]
  
#---------------------------generated covariates-------------------------------#
   se1 = 19 ; se2 = 94 #random seed
   set.seed(c(se1,se2), kind="Marsaglia-Multicarry") #To ensure repeatability of the experiment
   x1 <- runif(N, min = 0, max = 1)#covariate U(0,1)

   X <- matrix(c(rep(1,N), x1), ncol=2, byrow=F) #regressor matrix for the median submodel
   Z <- matrix(c(rep(1,N)), ncol = 1, byrow=F) #regressor matrix for the sigma submodel
   W <- matrix(c(rep(1,N)), ncol=1, byrow=F) #regressor matrix for the nu submodel

  
  #---------Initializing vectors to store the estimates and statistic z--------#
  thetahat <- matrix(NA, NREP, kk1 + kk2 + kk3)
  se_thetahat <- matrix(NA, NREP, kk1 + kk2 + kk3)
  z_thetahat <- matrix(NA, NREP, kk1 + kk2 + kk3)
  
  #contadores
  cont <- 0 # Counter for Monte Carlo replicates 
  n_descartadas <- 0 # Number of discarded replicates
  
  #
  while(cont < NREP)
  {
    cont <- cont + 1
    perc <- cont/NREP
    
    
    #--------------------To print the progress of the simulation-------------------#
    if(perc == 0.25 || perc == 0.5 || perc == 0.75 || perc ==1) cat("Perc. Replic. MC =",perc*100,"%","\n")
    
    #sample generation
    SampleG <- FunSample2(n = N, X = X, Z = Z, W = W,  theta = c(VALPHA, VGAMA, VXI), linkmu="log",
                         linksigma="log", linknu="log")
    
    #summary(y)
    y <- SampleG$y
    mu <- SampleG$mu
    sigma <- SampleG$sigma
    nu <- SampleG$nu
    
   #---- Starting values via RDTED ----#

     fitaux <-  suppressWarnings(tryCatch(gamlss(formula = y ~ SampleG$X[,2], 
                       sigma.formula = ~ 1, 
                       family = DTED(mu.link = "log", sigma.link="log"),
                       method = RS(), 
                       control = gamlss.control(n.cyc = 200, trace = F)), error=function(e) {e} )) 


     if(is.error(fitaux)){#starts else in case of error in the starting values
     muStart <- rep(median(y),length(y))
     sigmaStart <- rep(1.0,length(y))
     nuStart <- rep(1.0,length(y))
     }
     else
     {
     muStart <- exp(as.vector(X%*%as.numeric(fitaux$mu.coefficients)))
     sigmaStart <- rep(exp(as.numeric(fitaux$sigma.coefficients)), length(y))
     nuStart <- rep(1.0,length(y))
     }


    
    fit_result <- tryCatch({
      n_descartadas <- n_descartadas + 1
      
      fitMLE <- gamlss(formula = y ~ SampleG$X[,2], 
                       sigma.formula = ~ 1, 
                       nu.formula = ~ 1, 
                       family = DTWD(mu.link = "log",
                       sigma.link="log", nu.link="log"), method = RS(), 
                       control = gamlss.control(n.cyc = 200, trace = F),
                       mu.start = muStart, sigma.start = sigmaStart, nu.start = nuStart)
      
      #warning of not convergence
      if (fitMLE$converged && all(!is.na(sqrt(diag(vcov(fitMLE)))))) {
        list(success = TRUE, fit = fitMLE, se = sqrt(diag(vcov(fitMLE))), convergence = fitMLE$converged)
        
        
      }else{
        message("iteration error ", cont, ": algorithm did not converge.")
        list(success = FALSE, convergence = FALSE)
      }
      
      
    }, error = function(e) {
      
      #if error, print and return NULL
      message("iteration error ", cont, ": ", conditionMessage(e))
      list(success = FALSE, fit = FALSE, se = "FALSE", convergence = FALSE)
    })
    
    
    
    #-----------------------parameter estimates------------------------------#
    #
    if (fit_result$success & fit_result$convergence) {
      fitMLE <- fit_result$fit
      std_errors <- fit_result$se
      
      #estimate regression model coefficients
      thetahat[cont,] <- c(as.numeric(fitMLE$mu.coefficients), 
                           as.numeric(fitMLE$sigma.coefficients), 
                           as.numeric(fitMLE$nu.coefficients))
      
      se_thetahat[cont,] <- c(as.numeric(std_errors[1:ncol(X)]),
                              as.numeric(std_errors[(ncol(X) + 1):(ncol(X) + ncol(Z))]),
                              as.numeric(std_errors[(ncol(X) + ncol(Z) + 1):(ncol(X) + ncol(Z) + ncol(W))]))
      
      z_thetahat[cont, ] <- (thetahat[cont,] - VTHETA) / se_thetahat[cont,]
      
      
    } else {
      cont <- cont - 1  #discard and re-simulate this replication
      next
    }  
    
    
  }#closes MC replicates
  
  
  
  
  
  
  #-------------------------------RESULTS------------------------------------#
  
  
  
  #-----------------------mean of estimates-----------------------#
  M_thetahat <- apply(thetahat,2, function(x) mean(x, na.rm = TRUE))
  
  
  #---------------------median of estimates-----------------------#
  Med_thetahat <- apply(thetahat,2,function(x) median(x, na.rm = TRUE))
  
  #--------------------mean of s.e. estimates---------------------#
  M_se_thetahat <- apply(se_thetahat,2, function(x) mean(x, na.rm = TRUE))
  
  #--------empirical null levels of Wald-type tests---------------#
  LEVEL_MLE <- apply(abs(z_thetahat) > 1.959964, 2, function(x) mean(x, na.rm = TRUE))
  
  #----------------------- bias of estimates ---------------------#
  BIAS_thetahat <- (M_thetahat - VTHETA)
  
  #------------------- root mean squared error of estimates -------------------#
  RMSE_thetahat <- sqrt(apply(apply( t(thetahat) - VTHETA, 1, function(x) x^2), 2, mean) )
  
  
  #---------------------------------Outputs------------------------------------#
  #
  out_estimates_MLE <- cbind(M_thetahat, 
                             M_se_thetahat, 
                             LEVEL_MLE, 
                             Med_thetahat,
                             BIAS_thetahat,
                             RMSE_thetahat)
  
  cat(" ========================================================================================== \n")
  cat(" TRUE VALUES ","\n")
  cat(" Sample size ",N ,"\n")
  
  #
  out_VTHETA <- rbind(VTHETA)
  colnames(out_VTHETA) <- c(str_c("alpha", seq(1, (ncol(X)))),
                            str_c("gama", seq(1, (ncol(Z)))),
                            str_c("xi", seq(1, (ncol(W)))))
  
  
  
  rownames(out_VTHETA) <- c("")
  print(out_VTHETA)
  
  #
  Summ_mu <- round(c(min(mu, na.rm = T), 
                     max(mu, na.rm = T), 
                     median(mu, na.rm = T), 
                     min(sigma, na.rm = T), 
                     max(sigma, na.rm = T), 
                     median(sigma, na.rm = T), 
                     max(sigma)/min(sigma),
                     min(nu, na.rm = T), 
                     max(nu, na.rm = T), 
                     median(nu, na.rm = T),
                     n_desc = n_descartadas - NREP,
                     perc_desc = (n_descartadas - NREP)/NREP),3)
  
  #
  out_Smu <- rbind(Summ_mu)
  colnames(out_Smu)=c("Mu_min", 
                      "Mu_max", 
                      "Mu_median", 
                      "Sigma_min", 
                      "Sigma_max", 
                      "Sigma_median", 
                      "Het_Intensity",
                      "Nu_min", 
                      "Nu_max", 
                      "Nu_median",
                      "n_discard",
                      "perc_discard")
  
  rownames(out_Smu)=c("")
  print(out_Smu)
  
  #
  cat(" ========================================================================================== \n")
  cat(" ====================================== MLE estimates ===================================== \n")
  
  
  output_mle = rbind(out_estimates_MLE)
  colnames(output_mle)=c("Mean estimates","Mean SE", "Null level", "Median estimates","Bias estimates","RMSE estimates")
  rownames(output_mle)=c(str_c("MLE_alpha", seq(1, (ncol(X)))), 
                         str_c("MLE_gama", seq(1, (ncol(Z)))), 
                         str_c("MLE_xi", seq(1, (ncol(W)))))
  print(output_mle)
  
  cat(" ========================================================================================== \n")
  
  #------------------------------Print replicates--------------------------------#
  write.table(cbind(thetahat), 
              file = paste(names_results1[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
              col.names = c(str_c("MLE_alpha", seq(1, (ncol(X)))), 
                            str_c("MLE_gama", seq(1, (ncol(Z)))), 
                            str_c("MLE_xi", seq(1, (ncol(W))))))
  
  write.table(cbind(se_thetahat), 
              file= paste(names_results2[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
              col.names = colnames(out_VTHETA) <- c(str_c("se_MLE_alpha", seq(1, (ncol(X)))), 
                                                    str_c("se_MLE_gama", seq(1, (ncol(Z)))), 
                                                    str_c("se_MLE_xi", seq(1, (ncol(W))))))
  

  
  
}#closes samples loop

sink()
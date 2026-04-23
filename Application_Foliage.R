
if (!require("pacman")) install.packages("pacman")
pacman::p_load("ggplot2","qqboxplot","GLMsData", "gamlss", "gamlss.ggplots")

#*********** Foliage data application of DTED regression model *********#

source("DTED_GAMLSS.R")
source("DTWD_GAMLSS.R")

data(lime); attach(lime); names(lime)

# Descritive analysis

# Foliage versus Age

dataF <- data.frame(Age, Foliage)
ggplot(dataF, aes(x=Age, y=Foliage)) +
  geom_point(size=2, shape=19)+
  xlab("Age")+ylab("Foliage")+
  theme_bw()+
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
    )

#ggsave("SP_Age_Foliage.jpeg", width = 800, height = 500, units = 'px')

# Foliage versus Origin

dataF <- data.frame(Origin, Foliage)
ggplot(dataF, aes(x=Origin, y=Foliage)) + 
    geom_boxplot(fill="slateblue", alpha=0.5) +
    xlab("Origin")+ylab("Foliage")+
  theme_bw()+
  theme( 
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.x = element_text(
      size = 16,
      margin = margin(t = 20)  
    )
    )

ggsave("BP_Origin_Foliage.jpeg", width = 800, height = 500, units = 'px' )


#---- RDTED Fitted model ----#

# Covariates only in the median

fit_RDTED1 <- gamlss(formula = Foliage~ Age+Origin,
             sigma.formula = ~1,family = DTED(mu.link = "log",
             sigma.link="log"), method = RS())
summary(fit_RDTED1)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_RDTED1, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_RDTED <- resid(fit_RDTED1) #quantile residuals

dataR <-  data.frame(y = rq_RDTED)
  ggplot(dataR, aes(y = rq_RDTED)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
  )

#---- RDTWD Fitted model ----#

# Covariates only in the median

#Starting values from RDTED1 fit

muStart <-  fitted(fit_RDTED1, what= "mu")
sigmaStart <- fitted(fit_RDTED1, what= "sigma")
nuStart <- rep(1.0, length(Foliage))


fit_RDTWD1 <- gamlss(formula = Foliage~ Age+Origin,
              sigma.formula = ~1, nu.formula = ~1,
              family = DTWD(mu.link = "log",
              sigma.link="log", nu.link="log"), method = RS(),
              mu.start = muStart, sigma.start = sigmaStart,
              nu.start = nuStart)

summary(fit_RDTWD1)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_RDTWD1, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_RDTWD <- resid(fit_RDTWD1) #quantile residuals

dataR <-  data.frame(y = rq_RDTWD)
  ggplot(dataR, aes(y = rq_RDTWD)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
  )


#---- Gamma Fitted model ----#

# Covariates only in the mean


fit_Ga1 <- gamlss(formula = Foliage~ Age+Origin, 
sigma.formula = ~1,
family = GA(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_Ga1)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_Ga1, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_Ga1 <- resid(fit_Ga1) #quantile residuals

dataR <-  data.frame(y = rq_Ga1)
  ggplot(dataR, aes(y = rq_Ga1)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
    )

#---- Inverse Gaussian Fitted model ----#

# Covariates only in the mean

fit_IG1 <- gamlss(formula = Foliage~ Age+Origin, sigma.formula = ~1,
family = IG(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_IG1)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_IG1, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_IG1 <- resid(fit_IG1) #quantile residuals

dataR <-  data.frame(y = rq_IG1)
  ggplot(dataR, aes(y = rq_IG1)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
    )


#---- Weibull Fitted model ----#

# Covariates only in the mean

fit_W1 <- gamlss(formula = Foliage~ Age+Origin, sigma.formula = ~1,
family = WEI3(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_W1)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_W1, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_W1 <- resid(fit_W1) #quantile residuals

dataR <-  data.frame(y = rq_W1)
  ggplot(dataR, aes(y = rq_W1)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
    )


#---- Lognormal Fitted model ----#

# Covariates only in the median

fit_LN1 <- gamlss(formula = Foliage~ Age+Origin, sigma.formula = ~1,
family = LOGNO2(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_LN1)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_LN1, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_LN1 <- resid(fit_LN1) #quantile residuals

dataR <-  data.frame(y = rq_LN1)
  ggplot(dataR, aes(y = rq_LN1)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
    )

#---- RDTED Fitted model ----#

# Covariates only in mu and sigma

fit_RDTED2 <- gamlss(formula = Foliage ~ Age+Origin, sigma.formula = ~ Age-1,
              family = DTED(mu.link = "log", sigma.link="sqrt"), method = RS())
summary(fit_RDTED2)


#--- Diagnostics ---#

# Worm plot

resid_wp(fit_RDTED2, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_RDTED2 <- resid(fit_RDTED2) #quantile residuals

dataR <-  data.frame(y = rq_RDTED2)
  ggplot(dataR, aes(y = rq_RDTED2)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
  )


#---- RDTWD Fitted model ----#

# Covariates only in mu and sigma


#Starting values from RDTED2 fit

#muStart <-  fitted(fit_RDTED2, what= "mu")
#sigmaStart <- fitted(fit_RDTED2, what= "sigma")
#nuStart <- rep(1.0, length(Foliage))

# Or Starting values from RDTWD1 fit

muStart <-  fitted(fit_RDTWD1, what= "mu")
sigmaStart <- fitted(fit_RDTWD1, what= "sigma")
nuStart <- fitted(fit_RDTWD1, what= "nu")

fit_RDTWD2 <- gamlss(formula = Foliage ~ Age+Origin, sigma.formula = ~ Age-1,
              nu.formula=~1,family = DTWD(mu.link = "log", sigma.link="sqrt",
              nu.link="log"), method = RS(), mu.start = muStart,
              sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_RDTWD2)


#--- Diagnostics ---#

# Worm plot

resid_wp(fit_RDTWD2, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_RDTWD2 <- resid(fit_RDTWD2) #quantile residuals

dataR <-  data.frame(y = rq_RDTWD2)
  ggplot(dataR, aes(y = rq_RDTWD2)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
  )


#---- Gamma Fitted model ----#

# Covariates only in mu and sigma

fit_Ga2 <- gamlss(formula = Foliage~ Age+Origin, sigma.formula = ~Age-1,
family = GA(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_Ga2)


#--- Diagnostics ---#

# Worm plot

resid_wp(fit_Ga2, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_Ga2 <- resid(fit_Ga2) #quantile residuals

dataR <-  data.frame(y = rq_Ga2)
  ggplot(dataR, aes(y = rq_Ga2)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
    )

#---- Inverse Gaussian Fitted model ----#

# Covariates only in mu and sigma

fit_IG2 <- gamlss(formula = Foliage~ Age+Origin, sigma.formula = ~Age-1,
family = IG(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_IG2)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_IG2, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_IG2 <- resid(fit_IG2) #quantile residuals

dataR <-  data.frame(y = rq_IG2)
  ggplot(dataR, aes(y = rq_IG2)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
    )



#---- Weibull Fitted model ----#

# Covariates only in mu and sigma

fit_W2 <- gamlss(formula = Foliage~ Age+Origin, sigma.formula = ~Age-1,
family = WEI3(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_W2)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_W2, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_W2 <- resid(fit_W2) #quantile residuals

dataR <-  data.frame(y = rq_W2)
  ggplot(dataR, aes(y = rq_W2)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
    )


#----Lognormal Fitted model ----#

# Covariates only in mu and sigma

fit_LN2 <- gamlss(formula = Foliage~ Age+Origin, sigma.formula = ~Age-1,
family = LOGNO2(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_LN2)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_LN2, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_LN2 <- resid(fit_LN2) #quantile residuals

dataR <-  data.frame(y = rq_LN2)
  ggplot(dataR, aes(y = rq_LN2)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
    )


#---- RDTWD Fitted model ----#

# Covariates in mu, sigma and nu

#Starting values from RDTWD2 fit

muStart <-  fitted(fit_RDTWD2, what= "mu")
sigmaStart <- fitted(fit_RDTWD2, what= "sigma")
nuStart <- fitted(fit_RDTWD2, what= "nu")


fit_RDTWD3 <- gamlss(formula = Foliage ~ Age+Origin, sigma.formula = ~ Age-1,
              nu.formula=~Origin,family = DTWD(mu.link = "log", sigma.link="sqrt",
              nu.link="log"), method = RS(), mu.start = muStart,
             sigma.start = sigmaStart, nu.start = nuStart)
summary(fit_RDTWD3)


#--- Diagnostics ---#

# Worm plot

resid_wp(fit_RDTWD3, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_RDTWD3 <- resid(fit_RDTWD3) #quantile residuals

dataR <-  data.frame(y = rq_RDTWD3)
  ggplot(dataR, aes(y = rq_RDTWD3)) +                       
  geom_qqboxplot(notch=TRUE, varwidth = TRUE, reference_dist="norm") +
  xlab("reference: normal distribution") +
  ylab("Residuals") +
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(colour = "grey70"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    axis.text = element_text(size = 14)
  )



if (!require("pacman")) install.packages("pacman")
pacman::p_load("ggplot2","qqboxplot","GLMsData", "gamlss", "gamlss.ggplots")

#*********** Rent data application *********#

source("DTED_GAMLSS.R")
source("DTWD_GAMLSS.R")

data(rent99); attach(rent99); names(rent99)

# Descritive analysis

# rent versus area

dataR <- data.frame(area, rent)
ggplot(dataR, aes(x=area, y=rent)) +
  geom_point(size=2, shape=19)+
  xlab("Area")+ylab("Rent")+
  theme_bw()+
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
    )

# rent versus yearc

dataR <- data.frame(yearc, rent)
ggplot(dataR, aes(x=yearc, y=rent)) +
  geom_point(size=2, shape=19)+
  xlab("Yearc")+ylab("Rent")+
  theme_bw()+
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
    )

# rent versus location

dataR <- data.frame(location, rent)
ggplot(dataR, aes(x=location, y=rent)) + 
    geom_boxplot(fill="slateblue", alpha=0.5) +
    xlab("Location")+ylab("Rent")+
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

# rent versus bath

dataR <- data.frame(bath, rent)
ggplot(dataR, aes(x=bath, y=rent)) + 
    geom_boxplot(fill="slateblue", alpha=0.5) +
    xlab("Bath")+ylab("Rent")+
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

# rent versus kitchen

dataR <- data.frame(kitchen, rent)
ggplot(dataR, aes(x=kitchen, y=rent)) + 
    geom_boxplot(fill="slateblue", alpha=0.5) +
    xlab("Kitchen")+ylab("Rent")+
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


# Histogram of rent

df <- data.frame(x=rent)
 ggplot(df,aes(x=x))+
    geom_histogram(aes(y = ..density..),colour='white',fill='#696969')+
    labs(x='Rent',y='Density')+
    scale_color_manual(values=colors)+
    theme_bw()+
    theme(axis.title.y=element_text(colour='black',size=16),
          axis.title.x=element_text(colour='black',size=16),
          axis.text=element_text(colour='black',size=14),
          panel.border=element_blank(),
          axis.line=element_line(colour='black'))



#**********************************************

#---- RDTED Fitted model ----#

fit_RDTED <- gamlss(formula = rent ~ area+yearc+location+bath+kitchen, 
sigma.formula = ~area+location,
family = DTED(mu.link = "log", sigma.link="log"), method = RS(),
control = gamlss.control(n.cyc = 40, trace = T))
summary(fit_RDTED)


# Worm plot

resid_wp(fit_RDTED, ylim=2.0) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_RDTED <- resid(fit_RDTED) #quantile residuals

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

fit_RDTWD <- gamlss(formula =  rent ~ area+yearc+location+bath+kitchen,
              sigma.formula = ~ area+location, 
              nu.formula = ~ 1 ,
              family = DTWD(mu.link = "log",
              sigma.link="log", nu.link="log"), method = RS(), 
              control = gamlss.control(n.cyc = 40, trace = T))
summary(fit_RDTWD)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_RDTWD, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_RDTWD <- resid(fit_RDTWD) #quantile residuals

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

fit_Ga <- gamlss(formula = rent ~ area+yearc+location+bath+kitchen,
 sigma.formula = ~ area+location,
family = GA(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_Ga)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_Ga, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_Ga <- resid(fit_Ga) #quantile residuals

dataR <-  data.frame(y = rq_Ga)
  ggplot(dataR, aes(y = rq_Ga)) +                       
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

fit_IG <- gamlss(formula = rent ~ area+yearc+location+bath+kitchen,
 sigma.formula = ~area+location,
family = IG(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_IG)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_IG, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_IG <- resid(fit_IG) #quantile residuals

dataR <-  data.frame(y = rq_IG)
  ggplot(dataR, aes(y = rq_IG)) +                       
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

fit_W <- gamlss(formula = rent ~ area+yearc+location+bath+kitchen,
 sigma.formula = ~area+location,
family = WEI3(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_W)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_W, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_W <- resid(fit_W) #quantile residuals

dataR <-  data.frame(y = rq_W)
  ggplot(dataR, aes(y = rq_W)) +                       
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


fit_LN <- gamlss(formula = rent ~ area+yearc+location+bath+kitchen, 
sigma.formula = ~area+location,
family = LOGNO2(mu.link = "log", sigma.link="log"), method = RS())
summary(fit_LN)

#--- Diagnostics ---#

# Worm plot

resid_wp(fit_LN, ylim=2) +
  labs(title = NULL)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.background = element_rect(fill="white")
    )


# qqboxplot 

rq_LN <- resid(fit_LN) #quantile residuals

dataR <-  data.frame(y = rq_LN)
  ggplot(dataR, aes(y = rq_LN)) +                       
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


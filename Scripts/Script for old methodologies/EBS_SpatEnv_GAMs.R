library(tidyverse)
library(mgcv)
library(glmulti)


#Code written by Erin Fedewa

#Question: Are shifts in snow crab centers of distribution and area occupied related to bottom temperatures
  #and cold pool extent in the EBS? 
#Objective: Quantify changes in snow crab spatial indices over time in relation to env variables 

setwd('//Nmfs/akc-kod/Research/Spatial Opie MS/Datasets')


#################################### Add crab data ####################################################################
dat <- read.csv("Spatial_Env_TS.csv")
head(dat)
attach(dat)

#Questions: 
#How do GAM's deal with autocorrelated variables? 
#See Mike GitHub old code
#GAMM to account for spatial and temporal correlations more appropriate? 
#Use AICc for model with temp vrs cold pool extent? 


##################GAM for D95 immature female vrs cold pool extent and BT ##########################

#######D95 vrs Cold Pool Extent######## (Followed script from Mixed Effects Models and Extensions in Ecology)
plot(CP_EXTENT, D95_IMMATURE_FEMALE, type = "p")
#cross-validation used to estimate the optimal amt of smoothing
M3 <- gam(D95_IMMATURE_FEMALE ~ s(CP_EXTENT, fx = FALSE)) 
plot(M3, se = TRUE)
summary(M3) #Explained deviance is 38.1%, smoother for CP extent is significant 
anova(M3)

plot(CP_EXTENT, D95_IMMATURE_FEMALE, type = "p")
I1 <- order(CP_EXTENT)
lines(CP_EXTENT[I1], M3pred$fit[I1], lty=1)
lines(CP_EXTENT[I1], M3pred$fit[I1]+2*M3pred$se[I1],lty=2)
lines(CP_EXTENT[I1], M3pred$fit[I1]-2*M3pred$se[I1],lty=2)

#Validation Plots
gam.check(M3)
plot(fitted(M3), resid(M3)) #No underlying patterns in residuals
plot(CP_EXTENT, resid(M3))
hist(resid(M3))


#######D95 vrs Bottom Temp########
plot(AVG_BT, D95_IMMATURE_FEMALE, type = "p")
M4<- gam(D95_IMMATURE_FEMALE~ s(AVG_BT, fx = FALSE, k=-1, bs = "cr"))
plot(M4, se = TRUE)
summary(M4) #Explained deviance is 45.8%, smoother for temp is significant 
#Generalized cross-validation scores in output very large- 

#Validation Plots
gam.check(M4)
plot(fitted(M4), resid(M4)) #No underlying patterns in residuals
plot(CP_EXTENT, resid(M4))
hist(resid(M4))

#At least for D95 models, lm() seems the way to go, in initial plots, no nonlinear 
  #relationships that may warrant use of GAM's 




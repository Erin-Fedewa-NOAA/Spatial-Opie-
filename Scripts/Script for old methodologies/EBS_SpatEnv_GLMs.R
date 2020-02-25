library(tidyverse)
library(glmulti)
library(nlme)
library(jtools)
library(car)



#Code written by Erin Fedewa

#Question: Are shifts in snow crab centers of distribution and area occupied related to bottom temperatures
#and cold pool extent in the EBS across size/sex categories? 
#Objective: Quantify changes in snow crab spatial indices over time in relation to env variables 

setwd('//Nmfs/akc-kod/Research/Spatial Opie MS/Datasets')


#################################### Add crab data ####################################################################
dat <- read.csv("Spatial_Env_bySizeSex.csv")
head(dat)
#Create random size/sex term and center response variables to account for diff size/sex intercepts
dat %>%
  unite(SizeSex, c("SEX", "SIZE")) %>% 
  group_by(SizeSex) %>%
  mutate(center.D95 = center(D95)) %>%
  mutate(center.COD = center(LAT_COD)) ->dat #could also scale, but this is standardizing variance 
head(dat)


###################Data Exploration################################

#Data distribution
hist(dat$AVG_BT)
hist(dat$CP_EXTENT) #Right-skewed 
hist(dat$D95) #Bimodal b/w male and females 
hist(dat$scale.D95)
hist(dat$LAT_COD)

plot(AVG_BT ~ CP_EXTENT, data=dat)
plot(LAT_COD ~ D95, data=dat) #Correlation between response variables 
vif(lm(D95~CP_EXTENT + AVG_BT, data=dat)) #collinearity b/w covariates 

# look at the distribution of response variables 
dat %>%  #run on dat without combining sizesex
  ggplot()+
  geom_density(aes(x = D95))+
  facet_grid(cols = vars(SEX), rows = vars(SIZE))
  # No reason to think that D95 error is not normally distributed
dat %>%
  ggplot()+
  geom_density(aes(x = LAT_COD))+
  facet_grid(cols = vars(SEX), rows = vars(SIZE))
  #Skewed distributions for both male and female- specify error structure in GLM? 

# look for autocorrelation
dat %>%
  group_by(SEX, SIZE) %>%
  mutate(md95 = mean(D95)) %>%
  ggplot()+
  geom_line(aes(x = Year, y = D95))+
  geom_hline(aes(yintercept = md95), col = 2)+
  facet_grid(cols = vars(SEX), rows = vars(SIZE))

# 1st order AC should suffice for all size/sex categories 
dat %>%
  filter(SEX == "FEMALE",
         SIZE == "MATURE") -> tmp
acf(tmp$D95)

############# GLS with AR1 autocorrelation structure############################ 

######Scaled D95 vrs CP Extent/Bottom Temp Model 
  #Base model OLS vrs GLS with AR1 comparisons using maximum likelihood ratio test  
mod_0 <- gls(scale.D95 ~  CP_EXTENT, method="ML", data=dat) #standard OLS- base model 
  plot(mod_0) #distribution of error looks okay- no need to specify variance structure 
  resid<-resid(mod_0, type="normalized") #extract normalized residuals 
  acf(resid) #autocorrelation of residuals! Add AR1 correlation structure 
mod_1 <- gls(scale.D95 ~  CP_EXTENT, correlation = corAR1(), method="ML", data=dat)
  anova(mod_0, mod_1) #Correlation structure improves AIC scores   
#Check whether first order autocorrelation is sufficient 
mod1a<-update(mod_1, correlation=corARMA(p=2))
  anova(mod_1, mod1a) #AR1 looks good 

#Covariate comparisons using maximum likelihood ratio test 
mod_1 <- gls(scale.D95 ~  CP_EXTENT, correlation = corAR1(), method="ML", data=dat)
mod_2 <- gls(scale.D95 ~  AVG_BT, correlation = corAR1(), method="ML", data=dat)
  anova(mod_1, mod_2) 
  aicc(mod_1)
  aicc(mod_2) #Bottom temp appears to explain more variation in D95 than cold pool extent 

#Reapply final model with REML estimation (default method in nlme)
mod_final <- gls(scale.D95 ~ AVG_BT, correlation = corAR1(), data=dat)
  summary(mod_final)
  plot(mod_final)
  resid<-resid(mod_final, type="normalized") #extract normalized residuals 
  fit<-fitted(mod_final)
  plot(fit, resid)
  hist(resid) #Looks good
  acf(resid) #Also looks good 
  
  #Plots
  boxplot(predict(mod_final) ~ AVG_BT, data = dat)
  
  pred.mm <- ggpredict(mod_final, terms = c("AVG_BT"))  # this gives overall predictions for the model
  pred.mm
  
  # Plot the model predicted fit----Can't plot raw D95 data with scaled D95 model fit 
  ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +    # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = dat,                      # adding the raw data (scaled values)
               aes(x = AVG_BT, y = D95, colour = SizeSex)) + 
    labs(y="Snow Crab Areal Extent", x= "Cold Pool Extent") +
    scale_color_discrete(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males")) +
    theme_bw() +theme(legend.title=element_blank())
  
  #Plot 
  g1<-ggplot(dat, aes(AVG_BT, D95, col=SizeSex)) +
    geom_point() +
    geom_smooth(method="lm", se=F) + #raw data smoothed
   # geom_line(aes(y=predict(mod_final), group=SizeSex), size=1)+  #vrs fitted model predictions (same slope)
    labs(y="Snow Crab Areal Extent", x= "Cold Pool Extent") +
    scale_color_discrete(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males")) +
    theme_bw() +theme(legend.title=element_blank())
  g1
  
  

######Scaled COD vrs CP Extent/Bottom Temp Model 
  
  #Base model OLS vrs GLS with AR1 comparisons using maximum likelihood ratio test  
mod_0 <- gls(scale.COD ~  CP_EXTENT, method="ML", data=dat) #standard OLS- base model 
  plot(mod_0) #distribution of error looks okay- no need to specify variance structure 
  resid<-resid(mod_0, type="normalized") #extract normalized residuals 
  acf(resid) #autocorrelation of residuals! Add AR1 correlation structure 
mod_1 <- gls(scale.COD ~  CP_EXTENT, correlation = corAR1(), method="ML", data=dat)
  anova(mod_0, mod_1) #Correlation structure improves AIC scores   
  #Check whether first order autocorrelation is sufficient 
  mod1a<-update(mod_1, correlation=corARMA(p=2))
  anova(mod_1, mod1a) #AR1 still okay to use?? 
  
  #Covariate comparisons using maximum likelihood ratio test 
  mod_1 <- gls(scale.COD ~  CP_EXTENT, correlation = corAR1(), method="ML", data=dat)
  mod_2 <- gls(scale.COD ~  AVG_BT, correlation = corAR1(), method="ML", data=dat)
  anova(mod_1, mod_2) 
  aicc(mod_1)
  aicc(mod_2) #Best-fit model with cold pool fixed effect  
  
  #Reapply final model with REML estimation (default method in nlme)
  mod_final <- gls(scale.COD ~ CP_EXTENT, correlation = corAR1(), data=dat)
  summary(mod_final)
  plot(mod_final)
  resid<-resid(mod_final, type="normalized") #extract normalized residuals 
  fit<-fitted(mod_final)
  plot(fit, resid)
  hist(resid) #Looks good
  acf(resid) #Lag at 3, 4 and 5 years in residuals
  
  #Plots
  boxplot(predict(mod_final) ~ CP_EXTENT, data = dat)
  
  pred.mm <- ggpredict(mod_final, terms = c("CP_EXTENT"))  # this gives overall predictions for the model
  pred.mm
  
  # Plot the model predicted fit  ---Can't plot raw D95 data with scaled D95 model fit 
  ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +    # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = dat,                      # adding the raw data (scaled values)
               aes(x = AVG_BT, y = D95, colour = SizeSex)) + 
    labs(y="Snow Crab Areal Extent", x= "Cold Pool Extent") +
    scale_color_discrete(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males")) +
    theme_bw() +theme(legend.title=element_blank())
  
  #Plot raw data
  g1<-ggplot(dat, aes(CP_EXTENT, LAT_COD, col=SizeSex)) +
    geom_point() +
    geom_smooth(method="lm", se=F) + #raw data smoothed
    # geom_line(aes(y=predict(mod_final), group=SizeSex), size=1)+  #vrs fitted model predictions (same slope)
    labs(y="Snow Crab Center of Distribution", x= "Cold Pool Extent") +
    scale_color_discrete(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males")) +
    theme_bw() +theme(legend.title=element_blank())
  g1
  
  
  
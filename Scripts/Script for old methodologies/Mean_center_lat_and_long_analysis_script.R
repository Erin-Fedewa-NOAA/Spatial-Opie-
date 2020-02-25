##rm(list=ls())
################################## load libraries ##########################################################
library(ade4)
library(stats)
library(chron)
library(plyr)
library(ggplot2)
library(MASS)
library(car)
library(dplyr)
library(pander)
library(magrittr)
library(tables)
library(HSAUR)
library(graphics)
library(nlme)
library(Hmisc)
#####################################################################################################################################
################################## Correlations using Pyper and Peterman code #######################################################
#####################################################################################################################################


# Required function for cor.test.PP:

N.effective <- function(mat) {
# written by Franz Mueter   Last modified: 23 October 2000
# function to compute effective sample size for pairwise correlations among 
# autocorrelated variables. 
# Based on Pyper & Peterman (1998). CJFAS 55:2127-2140.  Eq. (1)
# where summation was done over j = 1, ..., N/5
# 
# mat a matrix of variables, one variable per column
#
# function to compute simple estimates of autocorrelation up to lag N/5: 
# (Eq. 7 in Pyper & Peterman)
 ar.fun <- function(x, max.lag = ceiling(sum(!is.na(x))/5)) {
   res <- rep(NA, max.lag)
   n <- length(x)
   for(i in 1.:max.lag) {
     x.bar <- mean(x, na.rm = T)
     res[i] <- ((n/(n - i)) * sum((x[1:(n - i)] - x.bar) * (x[(i + 1):n] - x.bar), na.rm
     = T))/sum((x - x.bar)^2, na.rm = T)
   }
   res
 }
 AR <- apply(mat, 2., ar.fun)
 k <- ncol(mat)
 if(is.matrix(AR)) {
   AR1 <- vector("list", k)
   for(i in 1:k) AR1[[i]] <- AR[, i]
   AR <- AR1  
 }
 N <- t(!is.na(mat)) %*% (!is.na(mat))
 N.lags <- ceiling(N/5.)
 N.eff <- matrix(0., k, k)
 # constrain effective N to smaller than or equal to actual N:
 for(i in 1.:k) {
   for(j in 1.:i) {
     lags <- 1.:N.lags[i, j]
     Nij <- N[i, j]
     N.eff[i, j] <- round((1./Nij + 2./Nij * sum((Nij - lags)/Nij * AR[[i]][lags] * AR[[
     j]][lags]))^(-1.))
   }
 }
 j <- N.eff > N
 N.eff[j] <- N[j]
 N.eff + t(N.eff) - diag(diag(N.eff))
}

# Function to test for significant correlation between two time series
# in the presence of autocorrelation in one or both of the series:

cor.test.PP <- function(x, y) {
# Function to test for significant correlations between x and y, which may be autocorrelated, 
# using modified Chelton method after Pyper & Peterman (1998)
# Eqn. 3 with N*-2 degrees of freedom
 N.eff <- N.effective(cbind(x, y))[1, 2]
 r <- cor(x, y, use="pair")
 fun <- function(alpha, N, r) {
   t2 <- qt(1 - alpha/2, N - 2)^2
   sqrt(t2/(t2 + N - 2)) - abs(r)
 }
 p.value <- uniroot(fun, c(1e-015, 0.9999), N = N.eff, r = r)$root
 cat("Two-sided test\n\n")
 c(correlation = r, P.value = p.value)
}



##################################################################################################################################################
################################### Add latitude and logitude for mean center data ###############################################################
##################################################################################################################################################

################################### Set work directory ###########################################################################################

setwd('//Nmfs/akc-kod/Research/Richar/Opilio_NBS_EBS/Common/Mean centers')

################################### Add data ######################################################################################################

dat<-read.csv("opilio_mean_center_lats_and_longs.csv")
dat

##################################################################################################################################################
################################# Add area of extent data for cold pool and opilio ###############################################################
##################################################################################################################################################

################################### Set work directory ###########################################################################################

setwd('//Nmfs/akc-kod/Research/Richar/Opilio_NBS_EBS/Common/Station_catch_pct')

dat2<-read.csv("Coldpool_ext_vs_opilio_extent.csv")
dat2
###################################################################################################################################################
################################# Analyze mean centers vs cold pool area ##########################################################################
###################################################################################################################################################
a<-dat$SURVEY_YEAR
x<-dat2$Coldpool_ext

###################################################################################################################################################
################################# MALES LE30 ######################################################################################################

################################# LATITUDE ##########################################################################################################
y<-dat$lat_malesLE30


plot(y~x,pch=16,ylab = "Mean center latitude",xlab = "Cold pool area (nmi^2)",main = "Males, LE30")
cor.test(x,y)            #no autocorrelation
cor.test.PP(x,y)         #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)
################################# LONGITUDE #######################################################################################################
y<-dat$long_malesLE30

plot(y~x,pch=16,ylab = "Mean center longitude",xlab = "Cold pool area (nmi^2)",main = "Males, LE30")
          
cor.test(x,y)            #no autocorrelation
cor.test.PP(x,y)         #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

###################################################################################################################################################
################################# MALES LE31to60 ######################################################################################################

################################# LATITUDE ##########################################################################################################
y<-dat$lat_males31to60


plot(y~x,pch=16,ylab = "Mean center latitude",xlab = "Cold pool area (nmi^2)",main = "Males, 31to60")

cor.test(x,y)            #no autocorrelation
cor.test.PP(x,y)         #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

################################# LONGITUDE #######################################################################################################
y<-dat$long_males31to60

plot(y~x,pch=16,ylab = "Mean center longitude",xlab = "Cold pool area (nmi^2)",main = "Males, 31to60")
          
cor.test(x,y)            #no autocorrelation
cor.test.PP(x,y)         #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

###################################################################################################################################################
################################# MALES LE61to90 ######################################################################################################

################################# LATITUDE ##########################################################################################################
y<-dat$lat_males61to90


plot(y~x,pch=16,ylab = "Mean center latitude",xlab = "Cold pool area (nmi^2)",main = "Males, 61to90")

cor.test(x,y)              #no autocorrelation
cor.test.PP(x,y)           #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

################################# LONGITUDE #######################################################################################################
y<-dat$long_males61to90

plot(y~x,pch=16,ylab = "Mean center longitude",xlab = "Cold pool area (nmi^2)",main = "Males, 61to90")
          
cor.test(x,y)      #no autocorrelation
cor.test.PP(x,y)   #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

###################################################################################################################################################
################################# MALES LE91to120 ######################################################################################################

################################# LATITUDE ##########################################################################################################
y<-dat$lat_males91to120


plot(y~x,pch=16,ylab = "Mean center latitude",xlab = "Cold pool area (nmi^2)",main = "Males, 91to120")


cor.test(x,y)           #no autocorrelation
cor.test.PP(x,y)        #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

################################# LONGITUDE #######################################################################################################
y<-dat$long_males91to120

plot(y~x,pch=16,ylab = "Mean center longitude",xlab = "Cold pool area (nmi^2)",main = "Males, 91to120")
          
cor.test(x,y)      #no autocorrelation
cor.test.PP(x,y)   #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

###################################################################################################################################################
################################# IMMATURE FEMALES ######################################################################################################

################################# LATITUDE ##########################################################################################################
y<-dat$lat_immatfem


plot(y~x,pch=16,ylab = "Mean center latitude",xlab = "Cold pool area (nmi^2)",main = "Immature females")

cor.test(x,y)           #no autocorrelation
cor.test.PP(x,y)        #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

################################# LONGITUDE #######################################################################################################
y<-dat$long_immatfem

plot(y~x,pch=16,ylab = "Mean center longitude",xlab = "Cold pool area (nmi^2)",main = "Immature females")
          
cor.test(x,y)      #no autocorrelation
cor.test.PP(x,y)   #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

###################################################################################################################################################
################################# MATURE FEMALES ######################################################################################################

################################# LATITUDE ##########################################################################################################
y<-dat$lat_matfem


plot(y~x,pch=16,ylab = "Mean center latitude",xlab = "Cold pool area (nmi^2)",main = "Mature females")

cor.test(x,y)           #no autocorrelation
cor.test.PP(x,y)        #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

################################# LONGITUDE #######################################################################################################
y<-dat$long_matfem

plot(y~x,pch=16,ylab = "Mean center longitude",xlab = "Cold pool area (nmi^2)",main = "Mature females")
          
cor.test(x,y)      #no autocorrelation
cor.test.PP(x,y)   #autocorrelation adjusted

fit<-lm(y~x)
summary(fit)

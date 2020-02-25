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
################################### Add extent data ###############################################################
##################################################################################################################################################

################################### Set work directory ###########################################################################################

setwd('//Nmfs/akc-kod/Research/Spatial Opie MS/Station_catch_pct/D95 EBS+NBS')

################################### Add data ######################################################################################################

dat<-read.csv("D95_TS.csv")
dat


###################################################################################################################################################
################################### Select predictor data columns #################################################################################################################

a<-dat$Year
x<-dat$Coldpool_extent                # This variable provides extent of cold pool for EBS + NBS
x2<-dat$Coldpool_NBS_extent           # This variable provides extent of cold pool for NBS
x3<-(x-x2)		        
x3   
                                 # This variable provides extent of cold pool for EBS                  
###################################################################################################################################################
################################### Opilio LE30 top 25% extent  #################################################################################################################

y<-dat$d95_males_le30




######################### EBS + NBS #########################################################
plot(y~x,pch=16, ylab="Opilio max density extent",xlab = "Cold pool extent", main = "Males LE30, EBS + NBS")
cor.test.PP(y,x)

fit<-lm(y~x)

summary(fit)

###################################################################################################################################################
################################### Opilio 31to60 top 25% extent  #################################################################################################################

y<-dat$d95_males_31to60


######################### EBS + NBS #########################################################
plot(y~x,pch=16, ylab="Opilio max density extent",xlab = "Cold pool extent", main = "Males 31 to 60, EBS + NBS")
cor.test.PP(y,x)

fit<-lm(y~x)

summary(fit)


###################################################################################################################################################
################################### Opilio 61to90 top 25% extent  #################################################################################################################

y<-dat$d95_males_61to90





######################### EBS + NBS #########################################################
plot(y~x,pch=16, ylab="Opilio max density extent",xlab = "Cold pool extent", main = "Males 61 to 90, EBS + NBS")
cor.test.PP(y,x)

fit<-lm(y~x)

summary(fit)



###################################################################################################################################################
################################### Opilio 91to120 top 25% extent  #################################################################################################################

y<-dat$d95_males_91to120





######################### EBS + NBS #########################################################
plot(y~x,pch=16, ylab="Opilio max density extent",xlab = "Cold pool extent", main = "Males 91 to 120, EBS + NBS")
cor.test.PP(y,x)

fit<-lm(y~x)

summary(fit)


###################################################################################################################################################
################################### Opilio immature female top 25% extent  #################################################################################################################

y<-dat$d95_immatfem





######################### EBS + NBS #########################################################
plot(y~x,pch=16, ylab="Opilio max density extent",xlab = "Cold pool extent", main = "Immature females, EBS + NBS")
cor.test.PP(y,x)

fit<-lm(y~x)

summary(fit)



###################################################################################################################################################
################################### Opilio mature female top 25% extent  #################################################################################################################

y<-dat$d95_matfem



######################### EBS + NBS #########################################################
plot(y~x,pch=16, ylab="Opilio max density extent",xlab = "Cold pool extent", main = "Mature females, EBS + NBS")
cor.test.PP(y,x)

fit<-lm(y~x)

summary(fit)








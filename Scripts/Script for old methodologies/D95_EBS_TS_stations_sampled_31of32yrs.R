library(tidyverse)
library(maps)
library(mapdata)
library(foreign)
library(nlme)

#######################################################################################################################
########################################### Pyper Peterman procedure for resolving autocorrelation in data ############

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

#############################################################################################################################
############################################### NBT ##########################################################################
##########Note to get both plots to plot together must run all code here 1st and then just the top portion for BB and minus BB (haven't looked into why)

mat<-matrix(1:2, 2, 1)
mat
layout(mat)
par(mar=c(2, 3, 0, 0.5), oma=c(0, 0, 0, 0), mgp=c(1.75,0.5,0))


setwd("//Nmfs/akc-kod/Survey/EBS Shelf/2019/Tech Memo/RawHaulAndCrab") #########UPDATE EACH YEAR
sta<-read.csv("stock_stations.csv")
attach(sta)
Allstations<-sta$Allstations
BBstations<-sta$BBRKC
Pribstations<-sta$PribRKC
minusBB<-Allstations[! Allstations %in% BBstations]
#################################################################################
setwd("//Nmfs/akc-kod/Survey/EBS Shelf/2019/Tech Memo/RawHaulAndCrab") ###########UPDATE EACH YEAR
dat<-read.csv("haul_newtimeseries.csv")
attach(dat)
dat2<-data.frame(SURVEY_YEAR, GIS_STATION, GEAR_TEMPERATURE, SURFACE_TEMPERATURE, PERFORMANCE, HAUL_TYPE)

ALLBS<-subset(dat2, HAUL_TYPE=="3" & PERFORMANCE>="0")

bot3<-data.frame(aggregate(GEAR_TEMPERATURE ~ SURVEY_YEAR, data=ALLBS, FUN=mean))
bot3


###############################################################################################################################
############################################### Subset bottom temperature for >1986 ###########################################

bot86<-subset(bot3,SURVEY_YEAR>=1988)
nbt86<-bot86$GEAR_TEMPERATURE

##############################################################################################################################
############################################### Cold pool data ###############################################################

CP<-subset(ALLBS,GEAR_TEMPERATURE <2.0)
CP
n<-nrow(CP)
count<-as.matrix(rep(1, times=n))
CP1<-cbind(CP,count)
colnames(CP1)

sy<-as.matrix(c(1975:2019))
n2<-nrow(sy)
out<-list()
for(i in 1:nrow(sy)){
     d<-subset(CP1,CP1$SURVEY_YEAR ==sy[i,])
      d<-as.data.frame(d)
      df<-sum(d$count)
#colnames(df)<-c("NumberStations")
		out[[i]]<-df
}

lf<-data.frame(matrix(unlist(out), nrow=1*n2,byrow=T))
sa<-as.matrix(rep(401,times = n2))
g<-cbind(sy,lf,sa) 
colnames(g)<-c("Year","NumberColdPoolStations","StationArea")
g

as.data.frame(g)
cpa<-as.matrix(g$NumberColdPoolStations*g$StationArea)
g1<-cbind(g,cpa)
colnames(g1)<-c("Year","NumberColdPoolStations","StationArea","ColdPoolArea")
g<-as.data.frame(g1)
g<-subset(g,Year>=1988)
cpa<-g$ColdPoolArea
as.matrix(cpa)

#################################### Set work directory for Crab Data ######################################################
setwd('//Nmfs/akc-kod/Research/Richar/Opilio_NBS_EBS/Common/Station_catch_pct/D95')


#################################### Add crab data ####################################################################
dat <- read.csv("co_cpue_ebs_data.csv")
range(dat$SURVEY_YEAR)

nrow(as.matrix(c(1988:2019)))
test<-subset(dat,SURVEY_YEAR==1988)
nrow(test)

# check
head(dat)
str(dat)

unique(dat$SURVEY_YEAR)
###########################################################################################################################
#################################### males LE30 ###########################################################################
###########################################################################################################################
# make a total CPUE colum
dat$total.cpue <- dat$MALE_LE30

# for now, restrict to 1988-2019
dat <- dat %>%
  filter(SURVEY_YEAR >= 1988)


# identify stations sampled in every year
ff <- function(x) sum(!is.na(x)) # create a function that counts the non-missing cases
samples <- tapply(dat$total.cpue, list(dat$SURVEY_YEAR, dat$GIS_STATION), ff)
samples

# change NA to 0
change <- is.na(samples)
samples[change] <- 0

# and look at total # of samples
total <- colSums(samples)
total

# how many stations are sampled each year?
sum(total>=31) # 369!

# restrict to these stations sampled each year
keepers <- total[total>=31]
dat <- dat %>%
  filter(GIS_STATION %in% names(keepers))

# and check 
check <- data.frame(one=unique(dat$GIS_STATION)[order(unique(dat$GIS_STATION))], two=names(keepers))
check # looks good!

# now plot the mean lat/long for these 305 stations sampled in every year
mean.lat <- tapply(dat$LATITUDE, as.character(dat$GIS_STATION), mean) 
# changing lat to a character allows us to drop missing values that are included when lat is a factor!
mean.long <- tapply(dat$LONGITUDE, as.character(dat$GIS_STATION), mean) 

plot(360+mean.long, mean.lat) # change longitude to ºE to plot
map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3")

# now...loop through and select the stations making up 95% of cpue in a year
years <- unique(dat$SURVEY_YEAR) # get an object of all the years in the data
years <- sort(years) # put them in order to keep things clean

output <- data.frame() # make an object to hold results

for(i in 1:length(years)){ # loop through each year
 # i <- 1
  temp <- dat %>%
    filter(SURVEY_YEAR==years[i]) # select the data 
  
  temp <- temp %>%
    arrange(desc(total.cpue)) # sort by cpue (large:small)
  
  temp$proportion.cpue <- temp$total.cpue / sum(temp$total.cpue) # calculate the proportion of total cpue for each station
  
  temp$cumulative.cpue <- cumsum(temp$proportion.cpue) # get the cumulative proportion for each station
  
  # and save results!
  output <- rbind(output, 
                  data.frame(year=years[i], 
                             station=temp$GIS_STATION,
                             d95 = temp$cumulative.cpue <= 0.95)) # T if in d95, F if not!
}

colnames(output)
d<-subset(output,d95=="TRUE")
d

t<-years
years<-as.matrix(years)
length(years)
#out <- data.frame()

out<-matrix(0,nrow=nrow(years), ncol = 1)

for(i in 1:nrow(years)){
	 d2<-subset(d, d$year==years[i,])
	n<-nrow(d2)
      a<-n*401
      out[i]<-a
 }


as.matrix(out)
ts<-cbind(years,out)
colnames(ts)<-c("Year","Area_D95")


###############################################################################################################################
################################################ Analyses ###################################################################
dev.new()
da<-cbind(ts,cpa)
da<-as.data.frame(da)

plot(da$Area_D95~da$cpa,pch=16, ylab = "Males LE30 D95 area",xlab = "Cold pool area", main = "Males LE30")

###########################################################################################################################
#################################### males 31to60 ###########################################################################
############################################################################################################################
dev.new()
# make a total CPUE colum
colnames(dat)
dat$total.cpue <- dat$MALE_31TO60

# for now, restrict to 1986-2019
dat <- dat %>%
  filter(SURVEY_YEAR >= 1988)
######################################


# identify stations sampled in every year
ff <- function(x) sum(!is.na(x)) # create a function that counts the non-missing cases
samples <- tapply(dat$total.cpue, list(dat$SURVEY_YEAR, dat$GIS_STATION), ff)
samples

# change NA to 0
change <- is.na(samples)
samples[change] <- 0

# and look at total # of samples
total <- colSums(samples)
total


# how many stations are sampled each year?
sum(total>=31) # 305!

# restrict to these stations sampled each year
keepers <- total[total>=31]
dat <- dat %>%
  filter(GIS_STATION %in% names(keepers))

# and check 
check <- data.frame(one=unique(dat$GIS_STATION)[order(unique(dat$GIS_STATION))], two=names(keepers))
check # looks good!

# now plot the mean lat/long for these 305 stations sampled in every year
mean.lat <- tapply(dat$LATITUDE, as.character(dat$GIS_STATION), mean) 
# changing lat to a character allows us to drop missing values that are included when lat is a factor!
mean.long <- tapply(dat$LONGITUDE, as.character(dat$GIS_STATION), mean) 

plot(360+mean.long, mean.lat) # change longitude to ºE to plot
map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3")

# now...loop through and select the stations making up 95% of cpue in a year
years <- unique(dat$SURVEY_YEAR) # get an object of all the years in the data
years <- sort(years) # put them in order to keep things clean

output <- data.frame() # make an object to hold results

for(i in 1:length(years)){ # loop through each year
 # i <- 1
  temp <- dat %>%
    filter(SURVEY_YEAR==years[i]) # select the data 
  
  temp <- temp %>%
    arrange(desc(total.cpue)) # sort by cpue (large:small)
  
  temp$proportion.cpue <- temp$total.cpue / sum(temp$total.cpue) # calculate the proportion of total cpue for each station
  
  temp$cumulative.cpue <- cumsum(temp$proportion.cpue) # get the cumulative proportion for each station
  
  # and save results!
  output <- rbind(output, 
                  data.frame(year=years[i], 
                             station=temp$GIS_STATION,
                             d95 = temp$cumulative.cpue <= 0.95)) # T if in d95, F if not!
}

colnames(output)
d<-subset(output,d95=="TRUE")
d

t<-years
years<-as.matrix(years)
length(years)
d

#out <- data.frame()
out<-matrix(0,nrow=nrow(years), ncol = 1)

for(i in 1:nrow(years)){
	 d2<-subset(d, d$year==years[i,])
	n<-nrow(d2)
      a<-n*401
      out[i]<-a
 }


as.matrix(out)
ts<-cbind(years,out)
colnames(ts)<-c("Year","Area_D95")

ts
###############################################################################################################################
################################################ Analyses ###################################################################

da<-cbind(ts,cpa)
da<-as.data.frame(da)
dev.new()
plot(da$Area_D95~da$cpa,pch=16, ylab = "Males 31to60 D95 area",xlab = "Cold pool area", main = "Males 31to60")

y<-da$Area_D95
x<-da$cpa

acf(y,main = "ACF plot for males 31-60mm")
fit<-lm(y~x)

summary(fit)
cor.test.PP(y,x)


fit2<-gls(y~x,correlation = corAR1())
summary(fit2)
###########################################################################################################################
#################################### males 61to90  ###########################################################################
###########################################################################################################################
dev.new()
# make a total CPUE colum
dat$total.cpue <- dat$MALE_61TO90

# for now, restrict to 1988-2019
dat <- dat %>%
  filter(SURVEY_YEAR >= 1988)


# identify stations sampled in every year
ff <- function(x) sum(!is.na(x)) # create a function that counts the non-missing cases
samples <- tapply(dat$total.cpue, list(dat$SURVEY_YEAR, dat$GIS_STATION), ff)
samples

# change NA to 0
change <- is.na(samples)
samples[change] <- 0

# and look at total # of samples
total <- colSums(samples)
total

# how many stations are sampled each year?
sum(total>=31) # 369!

# restrict to these stations sampled each year
keepers <- total[total>=31]
dat <- dat %>%
  filter(GIS_STATION %in% names(keepers))

# and check 
check <- data.frame(one=unique(dat$GIS_STATION)[order(unique(dat$GIS_STATION))], two=names(keepers))
check # looks good!

# now plot the mean lat/long for these 305 stations sampled in every year
mean.lat <- tapply(dat$LATITUDE, as.character(dat$GIS_STATION), mean) 
# changing lat to a character allows us to drop missing values that are included when lat is a factor!
mean.long <- tapply(dat$LONGITUDE, as.character(dat$GIS_STATION), mean) 

plot(360+mean.long, mean.lat) # change longitude to ºE to plot
map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3")

# now...loop through and select the stations making up 95% of cpue in a year
years <- unique(dat$SURVEY_YEAR) # get an object of all the years in the data
years <- sort(years) # put them in order to keep things clean

output <- data.frame() # make an object to hold results

for(i in 1:length(years)){ # loop through each year
 # i <- 1
  temp <- dat %>%
    filter(SURVEY_YEAR==years[i]) # select the data 
  
  temp <- temp %>%
    arrange(desc(total.cpue)) # sort by cpue (large:small)
  
  temp$proportion.cpue <- temp$total.cpue / sum(temp$total.cpue) # calculate the proportion of total cpue for each station
  
  temp$cumulative.cpue <- cumsum(temp$proportion.cpue) # get the cumulative proportion for each station
  
  # and save results!
  output <- rbind(output, 
                  data.frame(year=years[i], 
                             station=temp$GIS_STATION,
                             d95 = temp$cumulative.cpue <= 0.95)) # T if in d95, F if not!
}

colnames(output)
d<-subset(output,d95=="TRUE")
d

t<-years
years<-as.matrix(years)
length(years)
out<-matrix(0,nrow=nrow(years), ncol = 1)

for(i in 1:nrow(years)){
	 d2<-subset(d, d$year==years[i,])
	n<-nrow(d2)
      a<-n*401
      out[i]<-a
 }


as.matrix(out)
ts<-cbind(years,out)
colnames(ts)<-c("Year","Area_D95")


###############################################################################################################################
################################################ Analyses ###################################################################
dev.new()
da<-cbind(ts,cpa)
da<-as.data.frame(da)

plot(da$Area_D95~da$cpa,pch=16, ylab = "Males 61to90 D95 area",xlab = "Cold pool area", main = "Males 61to90")
y<-da$Area_D95
x<-da$cpa
acf(y,main = "ACF plot for males 61-90mm")

fit<-lm(y~x)

summary(fit)
cor.test.PP(y,x)

fit2<-gls(y~x,correlation = corAR1())
summary(fit2)
###########################################################################################################################
#################################### males 91to120  ###########################################################################
###########################################################################################################################
dev.new()
# make a total CPUE colum
dat$total.cpue <- dat$MALE_91TO120

# for now, restrict to 1988-2019
dat <- dat %>%
  filter(SURVEY_YEAR >= 1988)


# identify stations sampled in every year
ff <- function(x) sum(!is.na(x)) # create a function that counts the non-missing cases
samples <- tapply(dat$total.cpue, list(dat$SURVEY_YEAR, dat$GIS_STATION), ff)
samples

# change NA to 0
change <- is.na(samples)
samples[change] <- 0

# and look at total # of samples
total <- colSums(samples)
total

# how many stations are sampled each year?
sum(total>=31) # 369!

# restrict to these stations sampled each year
keepers <- total[total>=31]
dat <- dat %>%
  filter(GIS_STATION %in% names(keepers))

# and check 
check <- data.frame(one=unique(dat$GIS_STATION)[order(unique(dat$GIS_STATION))], two=names(keepers))
check # looks good!

# now plot the mean lat/long for these 305 stations sampled in every year
mean.lat <- tapply(dat$LATITUDE, as.character(dat$GIS_STATION), mean) 
# changing lat to a character allows us to drop missing values that are included when lat is a factor!
mean.long <- tapply(dat$LONGITUDE, as.character(dat$GIS_STATION), mean) 

plot(360+mean.long, mean.lat) # change longitude to ºE to plot
map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3")

# now...loop through and select the stations making up 95% of cpue in a year
years <- unique(dat$SURVEY_YEAR) # get an object of all the years in the data
years <- sort(years) # put them in order to keep things clean

output <- data.frame() # make an object to hold results

for(i in 1:length(years)){ # loop through each year
 # i <- 1
  temp <- dat %>%
    filter(SURVEY_YEAR==years[i]) # select the data 
  
  temp <- temp %>%
    arrange(desc(total.cpue)) # sort by cpue (large:small)
  
  temp$proportion.cpue <- temp$total.cpue / sum(temp$total.cpue) # calculate the proportion of total cpue for each station
  
  temp$cumulative.cpue <- cumsum(temp$proportion.cpue) # get the cumulative proportion for each station
  
  # and save results!
  output <- rbind(output, 
                  data.frame(year=years[i], 
                             station=temp$GIS_STATION,
                             d95 = temp$cumulative.cpue <= 0.95)) # T if in d95, F if not!
}

colnames(output)
d<-subset(output,d95=="TRUE")
d

t<-years
years<-as.matrix(years)
length(years)
out<-matrix(0,nrow=nrow(years), ncol = 1)

for(i in 1:nrow(years)){
	 d2<-subset(d, d$year==years[i,])
	n<-nrow(d2)
      a<-n*401
      out[i]<-a
 }


as.matrix(out)
ts<-cbind(years,out)
colnames(ts)<-c("Year","Area_D95")


###############################################################################################################################
################################################ Analyses ###################################################################
dev.new()
da<-cbind(ts,cpa)
da<-as.data.frame(da)

plot(da$Area_D95~da$cpa,pch=16, ylab = "Males 91to120 D95 area",xlab = "Cold pool area", main = "Males 91to120")
y<-da$Area_D95
x<-da$cpa
acf(y)

fit<-lm(y~x)

summary(fit)
cor.test.PP(y,x)

###########################################################################################################################
#################################### Immature females  ###########################################################################
###########################################################################################################################
dev.new()
# make a total CPUE colum
colnames(dat)
dat$total.cpue <- dat$IMMATURE_FEMALE

# for now, restrict to 1988-2019
dat <- dat %>%
  filter(SURVEY_YEAR >= 1988)


# identify stations sampled in every year
ff <- function(x) sum(!is.na(x)) # create a function that counts the non-missing cases
samples <- tapply(dat$total.cpue, list(dat$SURVEY_YEAR, dat$GIS_STATION), ff)
samples

# change NA to 0
change <- is.na(samples)
samples[change] <- 0

# and look at total # of samples
total <- colSums(samples)
total

# how many stations are sampled each year?
sum(total>=31) # 369!

# restrict to these stations sampled each year
keepers <- total[total>=31]
dat <- dat %>%
  filter(GIS_STATION %in% names(keepers))

# and check 
check <- data.frame(one=unique(dat$GIS_STATION)[order(unique(dat$GIS_STATION))], two=names(keepers))
check # looks good!

# now plot the mean lat/long for these 305 stations sampled in every year
mean.lat <- tapply(dat$LATITUDE, as.character(dat$GIS_STATION), mean) 
# changing lat to a character allows us to drop missing values that are included when lat is a factor!
mean.long <- tapply(dat$LONGITUDE, as.character(dat$GIS_STATION), mean) 

plot(360+mean.long, mean.lat) # change longitude to ºE to plot
map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3")

# now...loop through and select the stations making up 95% of cpue in a year
years <- unique(dat$SURVEY_YEAR) # get an object of all the years in the data
years <- sort(years) # put them in order to keep things clean

output <- data.frame() # make an object to hold results

for(i in 1:length(years)){ # loop through each year
 # i <- 1
  temp <- dat %>%
    filter(SURVEY_YEAR==years[i]) # select the data 
  
  temp <- temp %>%
    arrange(desc(total.cpue)) # sort by cpue (large:small)
  
  temp$proportion.cpue <- temp$total.cpue / sum(temp$total.cpue) # calculate the proportion of total cpue for each station
  
  temp$cumulative.cpue <- cumsum(temp$proportion.cpue) # get the cumulative proportion for each station
  
  # and save results!
  output <- rbind(output, 
                  data.frame(year=years[i], 
                             station=temp$GIS_STATION,
                             d95 = temp$cumulative.cpue <= 0.95)) # T if in d95, F if not!
}

colnames(output)
d<-subset(output,d95=="TRUE")
d

t<-years
years<-as.matrix(years)
length(years)
out<-matrix(0,nrow=nrow(years), ncol = 1)

for(i in 1:nrow(years)){
	 d2<-subset(d, d$year==years[i,])
	n<-nrow(d2)
      a<-n*401
      out[i]<-a
 }


as.matrix(out)
ts<-cbind(years,out)
colnames(ts)<-c("Year","Area_D95")


###############################################################################################################################
################################################ Analyses ###################################################################
dev.new()
da<-cbind(ts,cpa)
da<-as.data.frame(da)

plot(da$Area_D95~da$cpa,pch=16, ylab = "Immature female D95 area",xlab = "Cold pool area", main = "Immature females")
y<-da$Area_D95
x<-da$cpa
acf(y,main = "ACF plot for immature females")

fit<-lm(y~x)

summary(fit)
cor.test.PP(y,x)

fit2<-gls(y~x,correlation = corAR1())
summary(fit2)

###########################################################################################################################
#################################### Mature females  ###########################################################################
###########################################################################################################################
dev.new()
# make a total CPUE colum
colnames(dat)
dat$total.cpue <- dat$MATURE_FEMALE

# for now, restrict to 1988-2019
dat <- dat %>%
  filter(SURVEY_YEAR >= 1988)


# identify stations sampled in every year
ff <- function(x) sum(!is.na(x)) # create a function that counts the non-missing cases
samples <- tapply(dat$total.cpue, list(dat$SURVEY_YEAR, dat$GIS_STATION), ff)
samples

# change NA to 0
change <- is.na(samples)
samples[change] <- 0

# and look at total # of samples
total <- colSums(samples)
total

# how many stations are sampled each year?
sum(total>=31) # 369!

# restrict to these stations sampled each year
keepers <- total[total>=31]
dat <- dat %>%
  filter(GIS_STATION %in% names(keepers))

# and check 
check <- data.frame(one=unique(dat$GIS_STATION)[order(unique(dat$GIS_STATION))], two=names(keepers))
check # looks good!

# now plot the mean lat/long for these 305 stations sampled in every year
mean.lat <- tapply(dat$LATITUDE, as.character(dat$GIS_STATION), mean) 
# changing lat to a character allows us to drop missing values that are included when lat is a factor!
mean.long <- tapply(dat$LONGITUDE, as.character(dat$GIS_STATION), mean) 

plot(360+mean.long, mean.lat) # change longitude to ºE to plot
map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3")

# now...loop through and select the stations making up 95% of cpue in a year
years <- unique(dat$SURVEY_YEAR) # get an object of all the years in the data
years <- sort(years) # put them in order to keep things clean

output <- data.frame() # make an object to hold results

for(i in 1:length(years)){ # loop through each year
 # i <- 1
  temp <- dat %>%
    filter(SURVEY_YEAR==years[i]) # select the data 
  
  temp <- temp %>%
    arrange(desc(total.cpue)) # sort by cpue (large:small)
  
  temp$proportion.cpue <- temp$total.cpue / sum(temp$total.cpue) # calculate the proportion of total cpue for each station
  
  temp$cumulative.cpue <- cumsum(temp$proportion.cpue) # get the cumulative proportion for each station
  
  # and save results!
  output <- rbind(output, 
                  data.frame(year=years[i], 
                             station=temp$GIS_STATION,
                             d95 = temp$cumulative.cpue <= 0.95)) # T if in d95, F if not!
}

colnames(output)
d<-subset(output,d95=="TRUE")
d

t<-years
years<-as.matrix(years)
length(years)
out<-matrix(0,nrow=nrow(years), ncol = 1)

for(i in 1:nrow(years)){
	 d2<-subset(d, d$year==years[i,])
	n<-nrow(d2)
      a<-n*401
      out[i]<-a
 }


as.matrix(out)
ts<-cbind(years,out)
colnames(ts)<-c("Year","Area_D95")


###############################################################################################################################
################################################ Analyses ###################################################################
dev.new()
da<-cbind(ts,cpa)
da<-as.data.frame(da)

plot(da$Area_D95~da$cpa,pch=16, ylab = "Mature female D95 area",xlab = "Cold pool area", main = "Mature females")
y<-da$Area_D95
x<-da$cpa
acf(y,main = "ACF plot for mature females")

fit<-lm(y~x)

summary(fit)
cor.test.PP(y,x)

fit2<-gls(y~x,correlation = corAR1())
summary(fit2)

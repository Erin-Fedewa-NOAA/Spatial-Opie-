******#Male and Female C.o. NBS/EBS Size Frequency Plots for 2010/17/18/19

# Code written by Erin Fedewa/Tyler Jackson
#Hyp: Increase in large, mature snow crab in recent years in NBS as conditions
  #become "more like" EBS and crab potentially move north? 

# setwd
#setwd("//Nmfs/akc-kod/Research/Spatial Opie MS/Datasets")
setwd("C:/Users/erin.fedewa/Work/NBS/Distribution MS/Datasets")

library(tidyverse) 
library(plotrix)
library(cowplot)
library(WRS2)

# data ----
fem <- read_csv("CO_FEM_EBSNBS_1MM.csv") 
head(fem)
male <- read_csv("CO_MALE_EBSNBS_1MM.csv")
head(male)

# assign better heading names
names(fem) <- c("year", "region", "size", "Immature", "Mature", "Total")
names(male) <- c("year", "region", "size", "Mature", "Immature", "Total")
head(fem)

*************************#Females**************************************************

#NBS female abundance by maturity and year
  fem %>%
  filter(region =="NBS")%>%
  pivot_longer(cols = c(Immature, Mature)) %>%  ## reshape data so that you stack immature and mature crabs
  rename(class = name,
         abundance = value) %>%
  mutate(abundance_mil = abundance / 1000000, 
         abundance_mil = ifelse(class == "Immature", -1 * abundance_mil, abundance_mil)) %>%
  
  ggplot()+
  geom_bar(aes(y=abundance_mil, x=size, fill=class), stat="identity", width = 1, color = "grey40")+
  geom_hline(yintercept = 0)+
  xlim(0, 100)+
  scale_y_continuous(breaks = seq(-1000, 1000, 100))+
  labs(y = "Abundance (millions)", x = "Carapace width (mm)", fill = NULL)+
  facet_wrap(~year)+
  theme_bw()+
  theme(panel.grid = element_blank())
  
  
##Female abundance by maturity, region and year
fem %>%
  pivot_longer(cols = c(Mature, Immature)) %>%  ## reshape data so that you stack immature and mature crabs
  rename(class = name,
         abundance = value) %>%
  mutate(abundance_mil = abundance / 1000000, 
         abundance_mil = ifelse(region == "EBS", -1 * abundance_mil, abundance_mil)) %>%
  
  ggplot()+
  geom_bar(aes(y=abundance_mil, x=size, fill=class), stat="identity", position = "stack", width = 1, color = "grey40")+
  geom_hline(yintercept = 0)+
  facet_wrap(~year, scales='free')+
  theme(strip.text = element_text(face="bold", size=9))+
  scale_x_continuous(limits=c(10,85))+
  scale_y_continuous(limits=c(-500,500), breaks = seq(-1000, 1000, 100))+
  labs(y = "Abundance (millions)", x = "Carapace width (mm)", fill = NULL)+
  annotate("text", x = 80, y = 400, label = "NBS")+
  annotate("text", x = 80, y = -400, label = "EBS")+
  theme_bw()+
  theme(panel.grid = element_blank())

#Double checking data for very low abundances of imm females in 2019 
fem %>%
  filter(year == 2019,
         region =="EBS") %>%
  ggplot()+
  geom_bar(aes(x = size, y = Immature/1000000), stat = "identity")

##Proportion NBS Females by year
fem %>%
  filter(region == "NBS", size >= 30) %>%
  group_by(year) %>%
  mutate(Sum=sum(Total),
         prob = Total/Sum) ->dat2

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#by year
ggplot(dat2, aes(x=size, y=prob)) +
  geom_bar(stat ="identity", color = "black", fill = "grey60")+
  labs(y = "Proportion", x = "Carapace width (mm)") +
  scale_x_continuous(breaks = seq(0, 200, 10), limits = c(30, 80))+
  facet_wrap(~year, nrow = 4)+
  theme_bw()
h1<-5

#overlapping
g1<-ggplot(dat2, aes(x=size, y=prob, fill=factor(year))) +
  geom_density(stat = "identity", position ="identity", alpha=0.4)+
  labs(y = "Proportion", x = "Carapace width (mm)") +
  scale_x_continuous(breaks = seq(0, 200, 10), limits = c(30, 75))+
  scale_fill_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw()+theme(legend.title=element_blank())+theme(legend.position =c(.85,.85))
g1
#*****************Summary stats for female size frequency distributions 
  
  #Total Abundance/StDev of mature females by year/region
  fem %>%
  group_by(year, region) %>%
    summarize(Total=sum(Mature), SD=sd(Mature))
  
  #Total Abundance of all females  by year/region 
  fem %>%
    pivot_longer(cols = c(Mature, Immature)) %>%  
    rename(class = name, abundance = value) %>%
    group_by(year, region) %>%
    filter(size>=30) %>%
    summarize(Total=sum(abundance), Mean=weighted.mean(size, abundance, na.rm=TRUE))

  #Proportion of mature females in total female abundance 
  fem %>%
    pivot_longer(cols = c(Mature, Immature)) %>%  
    rename(class = name,
           abundance = value) %>%
    group_by(year, region) %>%
    filter(size>=30) %>%
    mutate(tot_abund = sum(abundance), 
           proportion = abundance / tot_abund * 100) %>%
    filter(class == "Mature") %>%
    summarize(Total=sum(proportion))
  
  #Calculating 90th percentile 
  fem %>%
    filter(region =="NBS") %>%
    filter(size>=30) %>%
    mutate(abundance_mil = Total / 1000000)->fem2
  
  #Find the 90th percentile for each year...this script is horrible- find a better way to do this in dplyr!
  NBS2010<-subset(fem2, year==2010)
  head(NBS2010)
  hist(NBS2010$abundance_mil)
  
  kurtosis(NBS2010$Total) #-1.17, thin tailed distribution- platykurtic
  skewness(NBS2010$Total)
  sizerep<-rep(NBS2010$size, times=NBS2010$abundance_mil)
  sizerep
  quantile(sizerep, c(.90)) #48mm
  NBS2010_90<-subset(NBS2010, size>=48) #Create dataframe with only 90th percentile
  
  
  NBS2017<-subset(fem2, year==2017)
  kurtosis(NBS2017$Total) #-0.10- thin tailed distribution- platykurtic
  skewness(NBS2017$Total)
  sizerep<-rep(NBS2017$size, times=NBS2017$abundance_mil)
  quantile(sizerep, c(.90)) #44mm
  NBS2017_90<-subset(NBS2017, size>=44) #Create dataframe with only 90th percentile 
  
  
  NBS2018<-subset(fem2, year==2018)
  kurtosis(NBS2018$Total) #1.39- fat tailed distribution- leptokurtic
  skewness(NBS2018$Total)
  sizerep<-rep(NBS2018$size, times=NBS2018$abundance_mil)
  quantile(sizerep, c(.90)) #46mm
  NBS2018_90<-subset(NBS2018, size>=46) #Create dataframe with only 90th percentile
  
  
  NBS2019<-subset(fem2, year==2019)
  kurtosis(NBS2019$Total) #-.64- thin tailed distribution- platykurtic
  skewness(NBS2019$Total)
  sizerep<-rep(NBS2019$size, times=NBS2019$abundance_mil)
  quantile(sizerep, c(.90)) #53mm
  NBS2019_90<-subset(NBS2019, size>=53) #Create dataframe with only 90th percentile 
  
  ###two sample K-S test to compare size distributions  
  ks.test(NBS2010$Total, NBS2019$Total) #Sig
  ks.test(NBS2010$Total, NBS2018$Total) #sig
  ks.test(NBS2010$Total, NBS2017$Total)
  ks.test(NBS2017$Total, NBS2018$Total)
  ks.test(NBS2017$Total, NBS2019$Total) #Sig
  ks.test(NBS2018$Total, NBS2019$Total)

  quantileTest(NBS2019$Total, NBS2010$Total, alternative="greater", target.quantile = 0.90)
  

  #Bonferroni p-value= 0.05/6 analyses = 0.008
  
  
*****************************#Males******************************************************
             
#NBS male abundance by maturity and year
    male %>%
    filter(region =="NBS")%>%
    pivot_longer(cols = c(Immature, Mature)) %>%  ## reshape data so that you stack immature and mature crabs
    rename(class = name,
           abundance = value) %>%
    mutate(abundance_mil = abundance / 1000000) %>%
    
    ggplot()+
    geom_bar(aes(y=abundance_mil, x=size, fill=class), stat="identity", position = "stack", width = 1, color = "grey40")+
    geom_hline(yintercept = 0)+
    facet_wrap(~year, scales='free')+
    scale_x_continuous(limits=c(0,100))+
    scale_y_continuous(limits=c(0,600))+
    labs(y = "Abundance (millions)", x = "Carapace width (mm)", fill = NULL)+
    geom_vline(xintercept = 78, linetype="longdash") +
        theme_bw()+
    theme(panel.grid = element_blank())
    
    
    
#Mature Male abundance by year and region
    male %>%
    mutate(abundance_mil = Mature / 1000000, 
           abundance_mil = ifelse(region == "EBS", -1 * abundance_mil, abundance_mil)) %>%
    
    ggplot()+
    geom_bar(aes(y=abundance_mil, x=size, fill=region), stat="identity", width = 1, color = "grey40")+
    geom_hline(yintercept = 0)+
    xlim(25, 160)+
    #scale_y_continuous(breaks = seq(-1000, 1000, 100))+
    labs(y = "Abundance (millions)", x = "Carapace width (mm)", fill = NULL)+
    #annotate("text", x = 0, y = 100, label = "NBS")+
    #annotate("text", x = 0, y = -100, label = "EBS")+
    geom_vline(xintercept = 78) +
    facet_wrap(~year, scales = "free_y" )+
    theme_bw()+
    theme(panel.grid = element_blank())  
    
    

    ##Male abundance by maturity, region and year 
  male %>%
  pivot_longer(cols = c(Mature, Immature)) %>%  ## reshape data so that you stack immature and mature crabs
  rename(class = name,
         abundance = value) %>%
  mutate(abundance_mil = abundance / 1000000, 
         abundance_mil = ifelse(region == "EBS", -1 * abundance_mil, abundance_mil)) %>%
  
  ggplot()+
  geom_bar(aes(y=abundance_mil, x=size, fill=class), stat="identity", position = "stack", width = 1, color = "grey40")+
  geom_hline(yintercept = 0)+
  facet_wrap(~year, scales='free')+
  theme(strip.text = element_text(face="bold", size=9))+
  scale_x_continuous(limits=c(10,120))+
  scale_y_continuous(limits=c(-300,600), breaks = seq(-1000, 1000, 100))+
  labs(y = "Abundance (millions)", x = "Carapace width (mm)", fill = NULL)+
  annotate("text", x = 115, y = 550, label = "NBS")+
  annotate("text", x = 115, y = -100, label = "EBS")+
  geom_vline(xintercept = 78, linetype="longdash") +
  theme_bw()+
  theme(panel.grid = element_blank())
  
  
#Proportion of total male abundance that are mature 
  male %>%
    pivot_longer(cols = c(Mature, Immature)) %>%  ## reshape data so that you stack immature and mature crabs
    rename(class = name,
           abundance = value) %>%
    group_by(year, region) %>%
        mutate(tot_abund = sum(abundance), 
           proportion = abundance / tot_abund * 100, 
           proportion = ifelse(region == "EBS", -1 * proportion, proportion)) %>%
       ungroup() %>%
    filter(class == "Mature") %>%
    
  ggplot()+
    geom_bar(aes(y=proportion, x=size, fill=region), stat="identity", width = 1, color = "grey40")+
    geom_hline(yintercept = 0)+
    xlim(0, 175)+
    labs(y = "Proportion of total male abundance (%)", x = "Carapace width (mm)", fill = NULL)+
    annotate("text", x = 0, y = 25, label = "NBS")+
    annotate("text", x = 0, y = -25, label = "EBS")+
    facet_wrap(~year, scales = "free_y" )+
    theme_bw()+
    theme(panel.grid = element_blank())
  
##Proportion NBS Males by year
  male %>%
    filter(region == "NBS", size>= 30) %>%
    group_by(year) %>%
    mutate(Sum=sum(Total),
           prob = Total/Sum) ->dat2

  #by year 
ggplot(dat2, aes(x=size, y=prob)) +
    geom_bar(stat ="identity", color = "black", fill = "grey60")+
    labs(y = "Proportion", x = "Carapace width (mm)") +
    scale_x_continuous(breaks = seq(0, 200, 10), limits = c(0, 100))+
    facet_wrap(~year, nrow = 4)+
    theme_classic()

#overlapping
g2<-ggplot(dat2, aes(x=size, y=prob, fill=factor(year))) +
  geom_density(stat = "identity", position ="identity", alpha=0.4, width=.85)+ #Kernel density estimates 
  labs(y = "Proportion", x = "Carapace width (mm)") +
  scale_x_continuous(breaks = seq(0, 200, 10), limits = c(30, 100))+
  scale_fill_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw()+theme(legend.title=element_blank())+theme(legend.position =c(.85,.85))
g2

#Male and Female NBS cowplot
plot_grid(g2, g1)#, labels = c('NBS Males', 'NBS Females'))
        
  
**************#Summary Stats********************************
#Total Abundance/StDev of mature males by year/region
male %>%
  group_by(year, region) %>%
  summarize(Total=sum(Mature), CI=(std.error(Mature)*1.96), Mean=weighted.mean(size, Mature, na.rm=TRUE))
  #Note: need by station abundance estimates for pooled CI's across mgmt area, not by size bin abundances 
            
 
#Total Abundance of all males by year/region 
  male %>%
    pivot_longer(cols = c(Mature, Immature)) %>%  
    rename(class = name, abundance = value) %>%
    group_by(year, region) %>%
    filter(size >= 30) %>%
    summarize(Total=sum(abundance), Mean=weighted.mean(size, abundance, na.rm=TRUE))
    

#Proportion of mature males in total male abundance 
male %>%
  pivot_longer(cols = c(Mature, Immature)) %>%  ## reshape data so that you stack immature and mature crabs
  rename(class = name,
         abundance = value) %>%
  group_by(year, region) %>%
  mutate(tot_abund = sum(abundance), 
         proportion = abundance / tot_abund * 100) %>%
  filter(class == "Mature") %>%
  summarize(Total=sum(proportion))

#Proportion of legal males in total male abundance 
male %>%
  pivot_longer(cols = c(Mature, Immature)) %>%  
  rename(class = name, 
         abundance = value) %>%
  group_by(year, region) %>%
  mutate(tot_abund = sum(abundance), 
         proportion = abundance / tot_abund * 100) %>%
  filter(size >= 78) %>%
  summarize(Total = sum(proportion))

#Calculating 90th percentile 
male %>%
  filter(region =="NBS") %>%
  filter(size>=30) %>%
  mutate(abundance_mil = Total / 1000000)->male2

#Find the 90th percentile for each year...this script is horrible- find a better way to do this in dplyr!
NBS2010<-subset(male2, year==2010)
kurtosis(NBS2010$Total) #7.96, fat tailed distribution- leptokurtic
sizerep<-rep(NBS2010$size, times=NBS2010$abundance_mil)
sizerep
quantile(sizerep, c(.90)) #51mm


NBS2017<-subset(male2, year==2017)
kurtosis(NBS2017$Total) #16.27- fat tailed distribution- leptokurtic
sizerep<-rep(NBS2017$size, times=NBS2017$abundance_mil)
quantile(sizerep, c(.90)) #45mm


NBS2018<-subset(male2, year==2018)
kurtosis(NBS2018$Total) #18.13- fat tailed distribution- leptokurtic
sizerep<-rep(NBS2018$size, times=NBS2018$abundance_mil)
quantile(sizerep, c(.90)) #47mm


NBS2019<-subset(male2, year==2019)
kurtosis(NBS2019$Total) #5.86- fat tailed distribution- leptokurtic
sizerep<-rep(NBS2019$size, times=NBS2019$abundance_mil)
quantile(sizerep, c(.90)) #66mm

###two sample K-S test to compare size distributions  
ks.test(NBS2010$Total, NBS2019$Total) #Sig
ks.test(NBS2010$Total, NBS2018$Total) 
ks.test(NBS2010$Total, NBS2017$Total)
ks.test(NBS2017$Total, NBS2018$Total)
ks.test(NBS2017$Total, NBS2019$Total) #Sig
ks.test(NBS2018$Total, NBS2019$Total)

#Bonferroni p-value= 0.05/6 analyses = 0.008



####################Combined Male/Female NBS Plots#####################################
# data ----
MF <- read_csv("CO_MALE&FEM_NBS_1MM.csv") 
head(MF)

##NBS abundance by sex, maturity and year
MF %>%
  pivot_longer(cols = c(ABUN_IMM, ABUN_MAT)) %>%  ## reshape data so that you stack immature and mature crabs
  rename(class = name, abundance = value) %>%
  mutate(abundance_mil = abundance / 1000000, 
         abundance_mil = ifelse(SEX == "FEMALE", -1 * abundance_mil, abundance_mil)) -> dat
  
  ggplot(dat)+
  geom_bar(aes(y=abundance_mil, x=SIZE1, fill=class), stat="identity", position = "stack", width = 1, color = "grey40")+
  geom_hline(yintercept = 0)+
  facet_wrap(~SURVEY_YEAR, scales='free_y')+
  theme(strip.text = element_text(face="bold", size=9))+
  scale_x_continuous(limits=c(10,100))+
  scale_y_continuous(breaks = seq(-1000, 1000, 100))+
  labs(y = "Abundance (millions)", x = "Carapace width (mm)", fill = NULL)+
  #annotate("text", x = 12, y = 450, label = "Male")+
  #annotate("text", x = 15, y = -450, label = "Female")+
  theme_bw()+
  theme(panel.grid = element_blank())


  
  

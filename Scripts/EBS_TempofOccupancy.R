# notes ----
#Calculate EBS snow crab temperatures of occupancy (CPUE weighted) from 1988-2019
#Cross correlations between temp of occupancy and bottom temps across time series 
#Rolling correlations to measure relationship between temp of occupancy/bottom temps over time, 
    #and to detent shifts in correlation over time
# Three general questions:
  # 1) how have avg. bottom temperatures changed in the EBS?
  # 2) how have bottom temperatures where oplio occur in the EBS changed?
  # 3) How has relationship b/w bottom temps and occupied temps changed over time? 

#Mike Litzow, Erin Fedewa, Tyler Jackson

# Last Updated 7/17/20

# load ----
library(tidyverse)
library(tidyquant)
library(stats)
library(TSA)
library(tseries)
library(forecast)
library(TTR)
library(FNGr) # for tickr function from Ben Williams, would need to download package from GitHub

# data ----

## EBS
catch_ebs <- read_csv("./Data/ebs_opilio_haul_data.csv")
head(catch_ebs)

# data mgmt ----

## crab data EBS
catch_ebs %>%
  filter(AKFIN_SURVEY_YEAR %in% c(1988:2019)) %>%
  select(AKFIN_SURVEY_YEAR, CRUISE, GIS_STATION, WIDTH, SHELL_CONDITION, SEX, CLUTCH_SIZE,
         SAMPLING_FACTOR, DISTANCE_FISHED, NET_WIDTH, MID_LATITUDE, MID_LONGITUDE,
         HAUL_TYPE, PERFORMANCE) -> crab_ebs
names(crab_ebs) <- c("year", "cruise", "Station", "cw", "sc", "sex", "clutch", 
                     "sample_factor", "distance_fished", "net_width", "lat", "lon", 
                     "haul_type", "performance")


## compute cpue by size-sex group for each station ----
crab_ebs %>% 
  filter(haul_type != 17,
         #performance == 0,
         cw >= 0) %>%
  mutate(size_sex = ifelse(sex == 1 & cw >= 31 & cw <= 60.9, "male31to60",
                           ifelse(sex == 1 & cw >= 61 & cw <= 90.9, "male61to90",
                                  ifelse(sex == 1 & cw >= 91 & cw <= 120.9, "male91to120",
                                         ifelse(sex == 2 & clutch >= 1, "mature_female",
                                                ifelse(sex == 2 & clutch == 0 & cw >= 31, "immature_female", NA))))),
         area_swept = distance_fished * (net_width / 1000) * 0.29155335) %>%
  group_by(year, Station, lat, lon, area_swept, size_sex) %>%
  summarise(num_crab = round(sum(sample_factor))) %>%
  filter(!is.na(area_swept)) %>%
  pivot_wider(names_from = size_sex, values_from = num_crab) %>%
  mutate(pop = sum(male31to60, male61to90, male91to120, immature_female, mature_female, na.rm = T)) %>%
  pivot_longer(c(6:12), names_to = "size_sex", values_to = "num_crab") %>%
  filter(size_sex != "NA") %>%
  mutate(num_crab = replace_na(num_crab, 0),
         cpue = num_crab / area_swept) %>%
  ungroup() -> cpue_long

#Merge imputed temps (from mice script) for temp occupancy calculations
dat2<-read.csv("./Data/Imputed Station Temps.csv")
head(dat2)

cpue_long %>%
  mutate(Station = gsub("-", "", Station)) %>%
  left_join(dat2, by = c("Station" = "GIS_STATION", "year" = "SURVEY_YEAR")) -> dat3

                
#Temperature of Occupancy Calculations ----
dat3 %>%
  filter(year >= 1988) %>%
  group_by(Station, size_sex) %>%
  filter(n() == 32) %>% #Include only stations sampled every year in 32 year timeseries 
  ungroup()%>% 
  group_by(year) %>%
  mutate(AVG_BT = mean(MEAN_GT)) %>%
  ungroup() %>%
  group_by(year, size_sex, AVG_BT) %>%
  summarise(TEMP_OCC = weighted.mean(MEAN_GT, w = cpue, na.rm = T)) -> dat3

#Output to include in master sizesex CSV 
write.csv(dat3, file="./Output/TempOcc_Output.csv")

# Cross-correlation analyses ----

# Functions to check for lags via cross-correlation
f_ccf <- function(x){ccf(x = x$AVG_BT, y = x$TEMP_OCC)}
f_max_acf <- function(x){max(x$acf)}
f_find_lag <- function(x){
  row <- which(x$acf == max(x$acf))
  x$lag[row]
}

#Max correlations/lags 
dat3  %>%
  nest(-size_sex) %>%
  mutate(ccf = purrr::map(data, f_ccf),
         max_cor = purrr::map_dbl(ccf, f_max_acf),
         lag = purrr::map_dbl(ccf, f_find_lag))

#Accounting for autocorrelation in ccf's ----

#Modified Chelton Method (run functions below first for cor.test.PP to work)
#Pyper and Peterman recommendation: use eq. 1 without the weighting
  #function, with autocorrelations estimated over N/5 lags j using
  #eq. 7, and with critical value (2) using N* - 2 degrees of freedom

dat3 %>%
  group_by(size_sex) %>%
  summarise(Corr = first(cor.test.PP(AVG_BT, TEMP_OCC)),
            P_val = last(cor.test.PP(AVG_BT, TEMP_OCC)))

########################################################

# Required functions for cor.test.PP: ----

N.effective <- function(mat) {
  # written by Franz Mueter   Last modified: 23 October 2000
  # function to compute effective sample size for pairwise correlations among
      # autocorrelated variables.
  # Based on Pyper & Peterman (1998). CJFAS 55:2127-2140.  Eq. (1)
     # where summation was done over lag j = 1, ..., N/5  #
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

################################################

#Plot 
dat3 %>%
  ggplot()+
  geom_line(aes(x = year, y = TEMP_OCC, col = size_sex), size=1)+
  geom_point(aes(x = year, y = TEMP_OCC,  col = size_sex), size=2)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y="Temperature (C)")+
  theme(axis.title.x = element_blank())+
  scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  theme(legend.title = element_blank())+
  theme(axis.text.x=element_text(size=11))+
  theme(axis.text.y=element_text(size=10)) +
  theme(legend.position = "bottom")
  #so there seems that there might be buffering during 2014-2017 in smaller size
    #classes, less so in 2018-2019....is this evident with rolling correlations? 

# 5-year rolling correlations b/w temps of occupancy and bottom temps ----
f_roll <- function(x){runCor(x = x$AVG_BT, y = x$TEMP_OCC, n=5)}

dat3 %>%
  group_by(size_sex)%>%
  nest() %>%
  mutate(roll = map(data, f_roll)) %>%
  unnest()  -> dat4

#Plot rolling correlations---can also present as 5 panel plot with static r for each

# set x axis labels
x_axis <- tickr(data = tibble(yr = 1988:2019), yr, 5)

# color palette
cbpalette <- c("#999999", "#F0E442", "#0072B2", "#009E73", "#D55E00", "black")

dat4 %>%
  filter(size_sex != "pop") %>%
  ggplot(aes(x = year, y = roll, group = size_sex, 
             shape = size_sex, color = size_sex))+
  geom_point()+
  geom_line()+
  scale_colour_manual(values = cbpalette, 
                      labels = c("Immature Females", "Mature Females", "Males 31 to 60 mm",
                                 "Males 61 to 90 mm", "Males 91 to 120 mm"))+
  scale_shape_manual(values = c(7:11),
                     labels = c("Immature Females", "Mature Females", "Males 31 to 60 mm",
                                "Males 61 to 90 mm", "Males 91 to 120 mm"))+
  labs(y = "Correlation Coefficient" , x = "",
       shape = NULL, color = NULL, alpha = NULL)+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(axis.text.x  = element_text(size=11)) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 9.5))
## write plot
ggsave(filename = "./Figs/5yr_corr.png", device = "png", width = 8, height = 6, 
       dpi = 300)
  
 #Linear models- are residuals autocorrelated, suggesting a temporal pattern not 
    #being captured by rxn b/w bottom temp and temp occ? 
dat3 %>%
  filter(size_sex == "immature_female") ->fem
 m1 <- lm(TEMP_OCC ~ AVG_BT, data=fem)
  summary(m1)
checkresiduals(m1)
#Repeated for all size.sex classes- seems more like white noise? 

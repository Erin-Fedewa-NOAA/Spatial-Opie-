
#notes ----
#Snow Crab Lat COD's in EBS 1988-2019 by size/sex category

#Erin Fedewa
#Last update 6/24/20

#load----
library(tidyverse)
library(rsample)

# data ----

#Exclude corner stations
corner <- list("QP2625","ON2625","HG2019","JI2120","IH1918",
               "GF2221","HG1918","GF2019","ON2524","PO2726",
               "IH2221","GF1918","JI2221","JI2019","JI1918",
               "HG2221","QP2726","PO2423","IH2019","PO2625",
               "QP2423","IH2120","PO2524","HG2120","GF2120",
               "QP2524")

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


## compute cpue by size-sex group for each station
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

#Compute EBS snow crab COD's by size/sex
cpue_long %>%
  filter(!(Station %in% corner)) %>% #exclude corner stations
  group_by(Station, size_sex) %>%
  filter(n() == 32) %>% #include only stations sampled every year in 32 year timeseries
  ungroup() %>% 
  group_by(year, size_sex) %>%
  summarise(Lat_COD = weighted.mean(lat, w = cpue)) -> COD 
    
write.csv(COD, file="./Output/COD_output.csv")
  








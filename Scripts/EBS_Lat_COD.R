#notes ----
#Lat/Long COD's in EBS 1988-2019
#.dbf files are output from GIS lat/long centroid calculations, weighted by CPUE 

#Erin Fedewa
#Last update 2/20/20

#load----
library(tidyverse)
library(foreign)

# data ----

#Lat/Long COD Timeseries ----

#ImmFem
dat <-read.dbf("./Data/Exported_shapefiles/EBS_immatfem_ge31.dbf")

dat %>%
  filter(SURVEY_YEA >= 1988) %>%
  select(SURVEY_YEA, LAT, LON) %>%  
  add_column(SizeSex = "FEMALE_IMM")->immfem

#MatFem
dat2 <-read.dbf("./Data/Exported_shapefiles/EBS_mature_females.dbf")

dat2 %>%
  filter(SURVEY_YEA >= 1988) %>%
  select(SURVEY_YEA, LAT, LONG) %>%  
  rename(LON = LONG) %>%
  add_column(SizeSex = "FEMALE_MAT")-> matfem

#Males 31-60
dat3 <-read.dbf("./Data/Exported_shapefiles/ebs_males_31to60s.dbf")

dat3 %>%
  filter(SURVEY_YEA >= 1988) %>%
  select(SURVEY_YEA, LAT, LONG) %>%  
  rename(LON = LONG) %>%
  add_column(SizeSex = "MALE_31TO60")-> male31_60

#Males 61-90
dat4 <-read.dbf("./Data/Exported_shapefiles/ebs_males_61to90s.dbf")

dat4 %>%
  filter(SURVEY_YEA >= 1988) %>%
  select(SURVEY_YEA, LAT, LONG) %>%  
  rename(LON = LONG) %>%
  add_column(SizeSex = "MALE_61TO90")-> male61_90

#Males 91 to 120
dat5 <-read.dbf("./Data/Exported_shapefiles/ebs_males_91to120s.dbf")

dat5 %>%
  filter(SURVEY_YEA >= 1988) %>%
  select(SURVEY_YEA, LAT, LONG) %>%  
  rename(LON = LONG) %>%
  add_column(SizeSex = "MALE_91TO120")-> male91_120

#Total Population
dat6 <-read.dbf("./Data/Exported_shapefiles/EBS_total_pop.dbf")

dat6 %>%
  filter(SURVEY_YEA >= 1988) %>%
  select(SURVEY_YEA, LAT, LON) %>%  
  add_column(SizeSex = "POP")-> pop

#Merge all dataframes and save output 
pop %>% 
  bind_rows(male31_60, male61_90, male91_120, immfem, matfem) -> COD

  write.csv(COD, file="./Output/Lat COD.csv") 
#Latitude COD's were imported into Spatial_Env_bySize_Sex.csv for LME/GLS models 

#Plot
ggplot(COD, aes(SURVEY_YEA, LAT, group=factor(SizeSex), colour=factor(SizeSex))) +
  geom_line() +
  facet_wrap(~SizeSex)

















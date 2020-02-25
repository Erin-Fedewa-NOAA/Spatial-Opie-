# notes ----
# Stratified Estimate of Abundance in NBS and EBS 1988 - 2019
#Fig 5 for Manuscript

# Tyler Jackson
# last updated: 2020/1/29

# load ----
library(tidyverse)
library(rsample)


# data ----

## NBS 
catch_nbs <- read_csv("./Data/nbs_opilio_haul_data.csv")
haul_nbs <- read_csv("./Data/haul_newtimeseries_nbs.csv")
strata_nbs <- read_csv("./Data/nbs_opilio_strata.csv")

## EBS
catch_ebs <- read_csv("./Data/ebs_opilio_haul_data.csv")
haul_ebs <- read_csv("./Data/haul_newtimeseries_ebs.csv")
strata_ebs <- read_csv("./Data/ebs_opilio_strata.csv")


# data mgmt ----

## crab data NBS
catch_nbs %>%
  full_join(haul_nbs, by = c("HAULJOIN", "GIS_STATION")) %>%
  select(SURVEY_YEAR, CRUISE.x, GIS_STATION, WIDTH, SHELL_CONDITION, SEX, CLUTCH_SIZE,
         SAMPLING_FACTOR, DISTANCE_FISHED, NET_WIDTH, MID_LATITUDE, MID_LONGITUDE,
         HAUL_TYPE, PERFORMANCE) -> crab_nbs
names(crab_nbs) <- c("year", "cruise", "Station", "cw", "sc", "sex", "clutch", 
                     "sample_factor", "distance_fished", "net_width", "lat", "lon", 
                     "haul_type", "performance")

## crab data EBS
catch_ebs %>%
  full_join(haul_ebs, by = c("HAULJOIN", "GIS_STATION")) %>%
  filter(AKFIN_SURVEY_YEAR %in% c(1988:2019)) %>%
  select(AKFIN_SURVEY_YEAR, CRUISE.x, GIS_STATION, WIDTH, SHELL_CONDITION, SEX, CLUTCH_SIZE,
         SAMPLING_FACTOR, DISTANCE_FISHED.x, NET_WIDTH.x, MID_LATITUDE.x, MID_LONGITUDE.x,
         HAUL_TYPE.x, PERFORMANCE.x) -> crab_ebs
names(crab_ebs) <- c("year", "cruise", "Station", "cw", "sc", "sex", "clutch", 
                     "sample_factor", "distance_fished", "net_width", "lat", "lon", 
                     "haul_type", "performance")

## combine crab data
crab <- bind_rows(crab_nbs, crab_ebs)

## strata NBS
names(strata_nbs) <- c("Station", "District", "tows", "stratum", "total_area", "year")

## strata EBS
strata_ebs %>%
  filter(SURVEY_YEAR %in% c(1988:2019)) %>%
  select(STATION_ID, DISTRICT, TOWS, STRATUM, TOTAL_AREA_SQ_NM, SURVEY_YEAR) -> strata_ebs
names(strata_ebs) <- c("Station", "District", "tows", "stratum", "total_area", "year")

## combine strata
strata <- bind_rows(strata_nbs, strata_ebs)

## compute cpue by size-sex group for each station
crab %>%
  filter(haul_type != 17,
         #performance == 0,
         cw >= 31) %>%
  mutate(size_sex = ifelse(sex == 1 & cw >= 31 & cw <= 60.9, "male31to60",
                           ifelse(sex == 1 & cw >= 61 & cw <= 90.9, "male61to90",
                                  ifelse(sex == 1 & cw >= 91 & cw <= 120.9, "male91to120",
                                         ifelse(sex == 2 & clutch >= 1, "mature_female",
                                                ifelse(sex == 2 & clutch == 0 & cw >= 31, "immature_female", NA))))),
         area_swept = distance_fished * (net_width / 1000) * 0.29155335)  %>%
  group_by(year, Station, lat, lon, area_swept, size_sex) %>%
  summarise(num_crab = round(sum(sample_factor)))  %>%
  right_join(strata, by = c("year", "Station")) %>%
  filter(!is.na(area_swept)) %>%
  pivot_wider(names_from = size_sex, values_from = num_crab) %>%
  pivot_longer(c(10:15), names_to = "size_sex", values_to = "num_crab") %>%
  filter(size_sex != "NA") %>%
  mutate(num_crab = replace_na(num_crab, 0),
         cpue = num_crab / area_swept,
         Region = ifelse(District == "NBS All", "NBS", "EBS")) %>%
  ungroup() -> cpue_long


haul_ebs %>%
  select(SURVEY_YEAR, GIS_STATION, STRATUM) %>%
  filter(SURVEY_YEAR %in% c(1988:2019)) %>%
  bind_rows(haul_ebs %>%
              select(SURVEY_YEAR, GIS_STATION, STRATUM) %>%
              filter(SURVEY_YEAR == 2019) %>%
              mutate(SURVEY_YEAR = 2018)) %>%
  bind_rows(haul_nbs %>%
              select(SURVEY_YEAR, GIS_STATION, STRATUM)) %>%
  rename(year = SURVEY_YEAR,
         Station = GIS_STATION,
         stratum_fish = STRATUM) %>%
  right_join(cpue_long, by = c("year", "Station")) %>%
  rename(stratum_crab = stratum) %>%
  mutate(stratum_crab = ifelse(stratum_crab == 999, 10, stratum_crab)) -> cpue_long


# get abundance by size sex catagory for 1988 - 2019 ----
cpue_long %>%
  group_by(year, Region, stratum_crab, size_sex) %>%
  summarise(total_area = mean(total_area),
            mean_cpue = mean(cpue),
            var_cpue = var(cpue),
            stratum_stations = n(),
            abundance_mil = mean(total_area) * mean_cpue / 1000000,
            var_abund_mil = mean(total_area)^2 * var_cpue / 1000000^2,
            stderr_mil = sqrt(var_abund_mil / n())) %>%
  group_by(year, Region, size_sex) %>%
  mutate(bs_area = sum(total_area),
         w = total_area / bs_area) %>% # define weights
  summarise(abundance_mil = sum(w * abundance_mil),
            stderr_mil = sum(w^2 * stderr_mil)) %>%
  #select(year, Region, size_sex, abundance_mil, stderr_mil) %>%
  mutate(lwr_95 = abundance_mil + qnorm(0.025) * stderr_mil,
         upp_95 = abundance_mil + qnorm(0.975) * stderr_mil) -> strat_est

## save EBS data as a .csv
strat_est %>%
  filter(Region != "NBS") %>%
  write_csv("./Output/strat_abund_ebs_1988_2019.csv")

##Save NBS data as a .csv
strat_est %>%
  filter(Region != "EBS") %>%
  write_csv("./Output/strat_abund_nbs.csv")

# THA plot! ----
## ebs
strat_est %>%
  filter(year %in% c(2010, 2017, 2018, 2019),
         Region == "EBS") %>%
  mutate(size_sex = case_when(size_sex == "male31to60" ~ "Males 31 - 60 mm",
                              size_sex == "male61to90" ~ "Males 61 - 90 mm",
                              size_sex == "male91to120" ~ "Males 91 - 120 mm",
                              size_sex == "immature_female" ~ "Immature Females",
                              size_sex == "mature_female" ~ "Mature Females"),
         size_sex = factor(size_sex, levels = c("Males 31 - 60 mm", "Males 61 - 90 mm", "Males 91 - 120 mm",
                                                "Immature Females", "Mature Females")),
         lwr_95 = ifelse(lwr_95 < 0, 0, lwr_95)) %>%
  ggplot()+
  geom_bar(aes(x = factor(year), y = abundance_mil), stat = "identity", fill = "navyblue", alpha = 0.7)+
  geom_errorbar(aes(x = factor(year), ymin = lwr_95, ymax = upp_95), width = 0.3)+
  geom_point(aes(x = factor(year), y = upp_95 * 1.1), color = NA)+
  labs(x = NULL, y = "Abundance (millions)", title = "EBS")+
  facet_wrap(~size_sex, scales = "free_y", ncol = 1)+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) -> ebs

## nbs
strat_est %>%
  filter(year %in% c(2010, 2017, 2018, 2019),
         Region == "NBS") %>%
  mutate(size_sex = case_when(size_sex == "male31to60" ~ "Males 31 - 60 mm",
                              size_sex == "male61to90" ~ "Males 61 - 90 mm",
                              size_sex == "male91to120" ~ "Males 91 - 120 mm",
                              size_sex == "immature_female" ~ "Immature Females",
                              size_sex == "mature_female" ~ "Mature Females"),
         size_sex = factor(size_sex, levels = c("Males 31 - 60 mm", "Males 61 - 90 mm", "Males 91 - 120 mm",
                                                "Immature Females", "Mature Females")),
         lwr_95 = ifelse(lwr_95 < 0, 0, lwr_95),) %>%
  ggplot()+
  geom_bar(aes(x = factor(year), y = abundance_mil), stat = "identity", fill = "navyblue", alpha = 0.7)+
  geom_errorbar(aes(x = factor(year), ymin = lwr_95, ymax = upp_95), width = 0.3)+
  geom_point(aes(x = factor(year), y = upp_95 * 1.1), color = NA)+
  labs(x = NULL, y = "Abundance (millions)", title = "NBS")+
  facet_wrap(~size_sex, scales = "free_y", ncol = 1)+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) -> nbs


ggsave("./Figs/abundance_bar_chart.png", plot = plot_grid(ebs, nbs), height = 6, width = 5, units = "in")







# notes ----
# K-S tests on NBS size frequency distributions by year
# Shift function to compare entire distributions (K-S test only testing central tendancy!)

# Tyler Jackson
# last updated: 2020/1/29

# load ----
library(tidyverse)

# data ----

## NBS 
catch_nbs <- read_csv("./Data/nbs_opilio_haul_data.csv")
haul_nbs <- read_csv("./Data/haul_newtimeseries_nbs.csv")
strata_nbs <- read_csv("./Data/nbs_opilio_strata.csv")


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

## strata NBS
names(strata_nbs) <- c("Station", "District", "tows", "stratum", "total_area", "year")

## compute cpue by size-sex group for each station
crab_nbs %>%
  mutate(cw_bin = round(cw, 0),
         area_swept = distance_fished * (net_width / 1000) * 0.539957^2) %>%
  group_by(year, Station, lat, lon, area_swept, sex, cw_bin) %>%
  summarise(num_crab = round(sum(sample_factor))) %>%
  ungroup() %>%
  right_join(strata_nbs, by = c("year", "Station")) %>%
  filter(!is.na(area_swept)) %>%
  unite(col = size_sex, sex:cw_bin, sep = "_") %>%
  pivot_wider(names_from = size_sex, values_from = num_crab) %>%
  pivot_longer(grep("\\d", names(.)), names_to = "size_sex", values_to = "num_crab") %>%
  select(-"NA_NA") %>%
  mutate(num_crab = replace_na(num_crab, 0),
         cpue = num_crab / area_swept) -> cpue_long

# stratified estimate of abundance in 1 mm bins ----
# NBS is only 1 stratied so there is no weighting. It is treated as a simple random sample

cpue_long %>%
  group_by(year, size_sex) %>%
  summarise(nbs_area = mean(total_area),
            mean_cpue = mean(cpue),
            var_cpue = var(cpue),
            abundance_mil = nbs_area * mean_cpue / 1000000,
            var_abund_mil = nbs_area^2 * var_cpue / 1000000^2,
            stderr_mil = sqrt(var_abund_mil / n()),
            lwr_95 = abundance_mil - 1.96 * stderr_mil,
            upp_95 = abundance_mil + 1.96 * stderr_mil) %>%
  select(year, size_sex, abundance_mil, lwr_95, upp_95) %>%
  separate(size_sex, sep = "_", into = c("sex", "size")) %>%
  arrange(size) %>%
  filter(size > 30) -> est_crab

## get "distribution" data
est_crab %>%
  select(1:4) %>%
  nest(-year, -sex) %>%
  mutate(dist = purrr::map(data, .f = function(x){as.numeric(with(x[x$size > 61,], 
                                                                  rep(size, abundance_mil)))})) -> tmp

# pairwise ks tests ----

## define a couple functions for ks test
## extract D statistic
f_stat <- function(x){x$statistic}
## extract p value
f_pval <- function(x){x$p.value}

### males   
tibble(y1 = c(2010, 2010, 2010, 2017, 2017, 2018),
       y2 = c(2017, 2018, 2019, 2018, 2019, 2019)) %>% #set up pairwise comparisons
  left_join(tmp %>%
              filter(sex == 1) %>%
              select(1, 4) %>%
              rename(dist1 = dist),
            by = c("y1" = "year")) %>%
  left_join(tmp %>%
              filter(sex == 1) %>%
              select(1, 4) %>%
              rename(dist2 = dist),
            by = c("y2" = "year")) %>% # add the distributions
  mutate(ks_test = purrr::map2(dist1, dist2, ks.test),
         D_stat = purrr::map_dbl(ks_test, f_stat),
         pval = purrr::map_dbl(ks_test, f_pval),
         pval_adj = p.adjust(pval, method = "bonferroni")) -> ks_male

#### to see the actual test output
ks_male %>%
  pull(ks_test) 

### females   
tibble(y1 = c(2010, 2010, 2010, 2017, 2017, 2018),
       y2 = c(2017, 2018, 2019, 2018, 2019, 2019)) %>% #set up pairwise comparisons
  left_join(tmp %>%
              filter(sex == 2) %>%
              select(1, 4) %>%
              rename(dist1 = dist),
            by = c("y1" = "year")) %>%
  left_join(tmp %>%
              filter(sex == 2) %>%
              select(1, 4) %>%
              rename(dist2 = dist),
            by = c("y2" = "year")) %>% # add the distributions
  mutate(ks_test = purrr::map2(dist1, dist2, ks.test),
         D_stat = purrr::map_dbl(ks_test, f_stat),
         pval = purrr::map_dbl(ks_test, f_pval),
         pval_adj = p.adjust(pval, method = "bonferroni")) -> ks_female

#### to see the actual test output
ks_female %>%
  pull(ks_test) 


# shift function ----
    # See https://garstats.wordpress.com/2016/07/12/shift-function/

## read shift function functions
source("./Scripts/Shift functions/Rallfun-v30.txt")
source("./Scripts/Shift functions/wilcox_modified.txt")
source("./Scripts/Shift functions/rgar_visualisation.txt")

### males 
est_crab %>%
  select(1:4) %>%
  nest(-year, -sex) %>%
  mutate(dist = purrr::map(data, 
                           .f = function(x){as.numeric(with(x[x$size > 61,], 
                                                            rep(size, abundance_mil)))})) -> tmp

tibble(y1 = c(2010, 2010, 2010, 2017, 2017, 2018),
       y2 = c(2017, 2018, 2019, 2018, 2019, 2019)) %>% #set up pairwise comparisons
  left_join(tmp %>%
              filter(sex == 1) %>%
              select(1, 4) %>%
              rename(dist1 = dist),
            by = c("y1" = "year")) %>%
  left_join(tmp %>%
              filter(sex == 1) %>%
              select(1, 4) %>%
              rename(dist2 = dist),
            by = c("y2" = "year")) %>% # add the distributions
  mutate(shift = purrr::map2(dist1, dist2, shifthd)) -> shift_m #name object at intermediate step since bootstrapping takes so long

shift_m %>%
  select(y1, y2, shift) %>%
  mutate(plot = purrr::map(shift, .f = function(x){
    x %>%
      ggplot()+
      geom_point(aes(x = group1, y = difference), col = "blue")+
      geom_line(aes(x = group1, y = difference), col = "blue")+
      geom_errorbar(aes(x = group1, ymin = ci_lower, ymax = ci_upper), width = 0, col = "blue")+
      geom_hline(yintercept = 0, linetype = 2)+
      labs(x = NULL, y = NULL)+
      #scale_x_continuous(limits = c(60, 75))+
      #scale_y_continuous(limits = c(-16, 6))+
      theme_bw()+
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
  })) -> x_m

### females

est_crab %>%
  select(1:4) %>%
  nest(-year, -sex) %>%
  mutate(dist = purrr::map(data, 
                           .f = function(x){as.numeric(with(x, 
                                                            rep(size, abundance_mil)))})) -> tmp


tibble(y1 = c(2010, 2010, 2010, 2017, 2017, 2018),
       y2 = c(2017, 2018, 2019, 2018, 2019, 2019)) %>% #set up pairwise comparisons
  left_join(tmp %>%
              filter(sex == 2) %>%
              select(1, 4) %>%
              rename(dist1 = dist),
            by = c("y1" = "year")) %>%
  left_join(tmp %>%
              filter(sex == 2) %>%
              select(1, 4) %>%
              rename(dist2 = dist),
            by = c("y2" = "year")) %>% # add the distributions
  mutate(shift = purrr::map2(dist1, dist2, shifthd)) -> shift_f

shift_f %>%
  select(y1, y2, shift) %>%
  mutate(plot = purrr::map(shift, .f = function(x){
    x %>%
      ggplot()+
      geom_point(aes(x = group1, y = difference), col = "red")+
      geom_line(aes(x = group1, y = difference), col = "red")+
      geom_errorbar(aes(x = group1, ymin = ci_lower, ymax = ci_upper), width = 0, col = "red")+
      geom_hline(yintercept = 0, linetype = 2)+
      labs(x = NULL, y = NULL)+
      #scale_x_continuous(limits = c(34, 50))+
      #scale_y_continuous(limits = c(-12, 6))+
      theme_bw()+
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
  })) -> x_f

cowplot::plot_grid(NULL, x_m$plot[[1]], x_m$plot[[2]], x_m$plot[[3]],
                   x_f$plot[[1]], NULL, x_m$plot[[4]], x_m$plot[[5]],
                   x_f$plot[[2]], x_f$plot[[4]], NULL,  x_m$plot[[6]],
                   x_f$plot[[3]], x_f$plot[[5]], x_f$plot[[6]], NULL,  
                   nrow = 4) -> x #add fig labels in pwpt

ggsave("./Figs/shift_example.png", plot = x, width = 8, height = 5, units = "in")

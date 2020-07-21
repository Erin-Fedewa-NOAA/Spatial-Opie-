# notes ----
# Pcod/Snow Crab Spatial Overlap in EBS and NBS 
# Erin Fedewa
# Last Updated 7/17/20

# load ----

library(tidyverse)
library(cowplot)
library(mgcv)
library(FNGr)
library(latticeExtra)

# data ----

#Add all EBS bottom trawl data files
ebs82 <- read_csv("./Data/Groundfish BT/ebs1982_1984.csv")
ebs85 <- read_csv("./Data/Groundfish BT/ebs1985_1989.csv")
ebs90 <- read_csv("./Data/Groundfish BT/ebs1990_1994.csv")
ebs95 <- read_csv("./Data/Groundfish BT/ebs1995_1999.csv")
ebs00 <- read_csv("./Data/Groundfish BT/ebs2000_2004.csv")
ebs05 <- read_csv("./Data/Groundfish BT/ebs2005_2008.csv")
ebs09 <- read_csv("./Data/Groundfish BT/ebs2009_2012.csv")
ebs13 <- read_csv("./Data/Groundfish BT/ebs2013_2016.csv")
ebs17 <- read_csv("./Data/Groundfish BT/ebs2017_2018.csv")
ebs19 <- read_csv("./Data/Groundfish BT/ebs2019.csv")

# combine datasets now and save output
bind_rows(ebs82, ebs85, ebs90, ebs95, ebs00, ebs05, ebs09, ebs13, ebs17, ebs19) %>%
  write_csv("./Output/ebs_timeseries.csv")

#Add NBS BT files 
nbs10 <- read_csv("./Data/Groundfish BT/nbs2010_17_19.csv")
nbs18 <- read_csv("./Data/Groundfish BT/nbs2018.csv") 

bind_rows(nbs10, nbs18) %>%
  write_csv("./Output/nbs_timeseries.csv")

# data managment ----

ebs <- read_csv("./Output/ebs_timeseries.csv")
head(ebs)
nbs <- read_csv("./Output/nbs_timeseries.csv")
head(nbs)

#Calculate % overlap for EBS timeseries 
ebs %>%
  group_by(YEAR) %>% 
  mutate(TOTAL_STATIONS = n_distinct(STATION)) %>%
  filter(SID %in% c("21720", "68580")) %>%
  select(YEAR, STATION, SID, WTCPUE, TOTAL_STATIONS) %>%
  mutate(SID = case_when(SID == 68580 ~ "CRAB",
                         SID == 21720 ~ "COD")) %>%
  group_by(YEAR, STATION, SID) %>%
  pivot_wider(names_from = SID, values_from = WTCPUE) %>%
  group_by(YEAR) %>%
  # method 1 -  % of total stations that include both cod and snow crab
  # method 2 - % of positive snow crab stations that included cod
  summarise(METHOD_1 = sum((CRAB > 0 & COD > 0), na.rm = T) / mean(TOTAL_STATIONS) * 100,
            METHOD_2 = sum((CRAB > 0 & COD > 0), na.rm = T) / sum((CRAB > 0), na.rm = T) * 100) -> overlap

 
#EBS Plot 
  # set x axis labels
  x_axis <- tickr(data = tibble(yr = 1980:2019), yr, 6)

  #Method 1
ggplot(aes(x = YEAR, y = METHOD_1), data = overlap) +
  geom_point() +
  geom_smooth(method = gam, formula = y~s(x, bs = "cs")) +
  labs(y = expression(atop("EBS Snow crab Pacific cod spatial overlap (%)")), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank())

#Method 2
ggplot(aes(x = YEAR, y = METHOD_2), data = overlap) +
  geom_point() +
  geom_smooth(method = gam, formula = y~s(x, bs = "cs")) +
  labs(y = expression(atop("EBS Snow crab Pacific cod spatial overlap (%)")), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank())
  

#Calculate % overlap for NBS timeseries
nbs %>%
  group_by(YEAR) %>% 
  mutate(TOTAL_STATIONS = n_distinct(STATIONID)) %>%
  filter(SPECIES_CODE %in% c("21720", "68580")) %>%
  select(YEAR, STATIONID, SPECIES_CODE, wCPUE, TOTAL_STATIONS) %>%
  mutate(SPECIES_CODE = case_when(SPECIES_CODE == 68580 ~ "CRAB",
                                  SPECIES_CODE == 21720 ~ "COD")) %>%
  group_by(YEAR, STATIONID, SPECIES_CODE) %>%
  pivot_wider(names_from = SPECIES_CODE, values_from = wCPUE) %>%
  group_by(YEAR) %>%
  # method 1 -  % of total stations that include both cod and snow crab
  # method 2 - % of positive snow crab stations that included cod
  summarise(METHOD_1 = sum((CRAB > 0 & COD > 0), na.rm = T) / mean(TOTAL_STATIONS) * 100,
            METHOD_2 = sum((CRAB > 0 & COD > 0), na.rm = T) / sum((CRAB > 0), na.rm = T) * 100) %>%
  add_row(YEAR = 2011:2016, METHOD_1 = NA, METHOD_2 = NA) -> nbs_overlap  #Add NA's for missing yrs 

#NBS Plot 
# set x axis labels
x_axis <- tickr(data = tibble(yr = 2008:2019), yr, 2)

  #Method 1
ggplot(aes(x = YEAR, y = METHOD_1), data = nbs_overlap, na.rm = T) +
  geom_point() +
  geom_line() +
  labs(y = expression(atop("NBS Snow crab Pacific cod spatial overlap (%)")), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank())

#Method 2
  ggplot(aes(x = YEAR, y = METHOD_2), data = nbs_overlap, na.rm = T) +
  geom_point(color = "#E69F00", size = 6) +
  geom_line() +
  labs(y = expression(atop("Snow crab-Pacific cod spatial overlap (%)")), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey85", size = 0.2)) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11)) -> overlap

# Load northern Bering Sea haul data
nbs <- read_csv("./Data/haul_newtimeseries_nbs.csv")

# compute nbs average temperatures
nbs %>%
  rename(YEAR = SURVEY_YEAR) %>%
  group_by(YEAR) %>%
  summarise(AVG_BT = mean(GEAR_TEMPERATURE)) %>%
  add_row(YEAR = 2011:2016, AVG_BT = NA) -> nbs_bt

#NBS temp Plot 
nbs_bt %>%
  select(YEAR, AVG_BT) %>%
  ggplot(aes(x = YEAR, y = AVG_BT, group = 1), na.rm = T)+
  geom_point(color = "#E69F00", size = 6)+
  geom_line()+
  #geom_hline(aes(yintercept = mean(AVG_BT)), linetype = 2) +
  labs(y = expression(atop("Mean Bottom Temperature "( degree~C))), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey85", size = 0.2)) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11)) -> avg_bt

#Combine NBS plots 
plot_grid(avg_bt, overlap + 
          theme(legend.position = "none"), 
          labels = c('(a)', '(b)'),
          label_size = 11,
          hjust = -5, vjust = 2.5,
          nrow = 2, align = "hv", axis = "l", rel_heights = c(1, 1)) -> pred_plot


#Double Y axis plot in ggplot ---NOT WORKING
nbs <- inner_join(nbs_overlap, nbs_bt, by = "YEAR") 
nbs %>%
  ggplot(aes(x=YEAR), na.rm=T) + 
  geom_point(aes(y=AVG_BT, colour = "Mean Temperature")) +
  geom_line(aes(y=AVG_BT, colour = "Spatial Overlap"))+
  #Add predator overlap data, transformed to match the range of temperature 
  geom_point(aes(y=METHOD_2/20, color = "#0072B2")) +
  geom_line(aes(y=METHOD_2/20, color = "#0072B2"))+
  scale_y_continuous(sec.axis = sec_axis(~.*20, name = "Snow crab-Pacific cod spatial overlap (%)"))

#Double Y axis plot in lattice 
Temp <- xyplot(AVG_BT ~ YEAR, nbs_bt, type = "o" , lwd=2, col="#E69F00", 
               ylab="Mean Bottom Temperature", xlab="",)
Pred <- xyplot(METHOD_2 ~ YEAR, nbs_overlap, type = "l", lwd=2, col="#0072B2", grid = TRUE,
               ylab="Snow crab-Pacific cod spatial overlap (%)", xlab="")
#Make the plot with second y axis:
  doubleYScale(Temp, Pred, add.ylab2 = TRUE, use.style=FALSE)


## write plot
ggsave(filename = "./Figs/pred_overlap.png", device = "png", width = 5, height = 7, 
       dpi = 300) 

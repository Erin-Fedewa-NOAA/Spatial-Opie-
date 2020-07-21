# notes ----
# Time Series Figure 2 in Manuscript 
#Temperature difference plot for Reviewer response
# Last Updated 7/17/20

# load ----

library(tidyverse)
library(cowplot)
library(FNGr) # for tickr function from Ben Williams, would need to download package from GitHub

# set x axis labels
x_axis <- tickr(data = tibble(yr = 1988:2019), yr, 5)

# color palette
cbpalette <- c("#999999", "#F0E442", "#0072B2", "#009E73", "#D55E00", "black")
# data ----

# all time series
ts <- read_csv("./Data/Spatial_Env_bySizeSex.csv")
head(ts)


# northern Bering Sea haul data
nbs <- read_csv("./Data/haul_newtimeseries_nbs.csv")

# compute nbs average temperatures
nbs %>%
  group_by(SURVEY_YEAR) %>%
  summarise(AVG_BT = mean(GEAR_TEMPERATURE)) -> nbs_bt

# Fig 2 for manuscript ----

## create the plot of environmental timeseries
ts %>%
  ## mean bottom temperature
  select(Year, AVG_BT) %>%
  ggplot(aes(x = Year, y = AVG_BT, group = 1))+
  geom_point()+
  geom_point(data = nbs_bt, aes(x = SURVEY_YEAR, y = AVG_BT), shape = 8, size = 1.5)+
  geom_line()+
  geom_hline(aes(yintercept = mean(AVG_BT)), linetype = 2)+
  labs(y = expression(atop("Mean Bottom", "Temperature "( degree~C))), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank()) -> avg_bt

ts %>%
  ## cold pool areal extent
  select(Year, CP_EXTENT) %>%
  ggplot(aes(x = Year, y = CP_EXTENT, group = 1))+
  geom_point()+
  geom_line()+
  geom_hline(aes(yintercept = mean(CP_EXTENT)), linetype = 2)+
  labs(y = bquote('Cold Pool Spatial Extent ('~nmi~')'), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank()) -> cpa

ts %>%
  ## cold pool center of distribution
  select(Year, CP_COD) %>%
  ggplot(aes(x = Year, y = CP_COD, group = 1))+
  geom_point()+
  geom_line()+
  geom_hline(aes(yintercept = mean(CP_COD)), linetype = 2)+
  labs(y = expression(atop("Cold Pool Center of", "Distribution "( degree~Latitude))), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank()) -> cp_cod


## create the plot of biological timeseries
ts %>%
  # Temperature of Occupancy
  select(Year, SIZESEX, TEMP_OCC)  %>%
  mutate(mean_Pop = mean(TEMP_OCC[.$SIZESEX == "POP"])) %>%
  ggplot(aes(x = Year, y = TEMP_OCC, group = SIZESEX, 
             alpha = SIZESEX, shape = SIZESEX, color = SIZESEX))+
  geom_point()+
  geom_line()+
  geom_hline(aes(yintercept = mean_Pop), linetype = 2)+
  #scale_alpha_manual(values = c(1, 1, 1, 1, 1, 1))+
  scale_alpha_manual(values = c(0.7, 0.7, 0.7, 0.7, 0.7, 1))+
  scale_colour_manual(values = cbpalette)+
  scale_shape_manual(values = c(7:11, 16))+
  labs(y = expression(atop("Temperature of", "Occupancy "( degree~C))), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none") -> temp_occ

ts %>%
  # D95
  select(Year, SIZESEX, D95)  %>%
  mutate(mean_Pop = mean(D95[.$SIZESEX == "POP"])) %>%
  ggplot(aes(x = Year, y = D95, group = SIZESEX, 
             alpha = SIZESEX, shape = SIZESEX, color = SIZESEX))+
  geom_point()+
  geom_line()+
  geom_hline(aes(yintercept = mean_Pop), linetype = 2)+
  scale_alpha_manual(values = c(0.7, 0.7, 0.7, 0.7, 0.7, 1))+
  scale_colour_manual(values = cbpalette)+
  scale_shape_manual(values = c(7:11, 16))+
  labs(y = bquote('Snow Crab Spatial Extent ('~nmi~')'), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none") -> d95

ts %>%
  # COD 
  select(Year, SIZESEX, LAT_COD) %>%
  mutate(mean_Pop = mean(LAT_COD[.$SIZESEX == "POP"])) %>%
  ggplot(aes(x = Year, y = LAT_COD, group = SIZESEX, 
             alpha = SIZESEX, shape = SIZESEX, color = SIZESEX))+
  geom_point()+
  geom_line()+
  geom_hline(aes(yintercept = mean_Pop), linetype = 2)+
  scale_alpha_manual(values = c(0.7, 0.7, 0.7, 0.7, 0.7, 1),
                     labels = c("Immature Females", "Mature Females", "Males 31 to 60 mm",
                                "Males 61 to 90 mm", "Males 91 to 120 mm", "Population"))+
  scale_colour_manual(values = cbpalette, 
                      labels = c("Immature Females", "Mature Females", "Males 31 to 60 mm",
                                 "Males 61 to 90 mm", "Males 91 to 120 mm", "Population"))+
  scale_shape_manual(values = c(7:11, 16),
                     labels = c("Immature Females", "Mature Females", "Males 31 to 60 mm",
                                "Males 61 to 90 mm", "Males 91 to 120 mm", "Population"))+
  labs(y = expression(atop("Snow Crab Center of", "Distribution "( degree~Latitude))) , x = "",
       shape = NULL, color = NULL, alpha = NULL)+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10)) -> cod

# extract legend 
key <- get_legend(cod)

# create an empty plot
ggplot()+
  theme(panel.background = element_blank(),
        panel.border = element_blank()) -> empty

## combine plots
plot_grid(avg_bt, temp_occ, 
          cpa, d95, 
          cp_cod, cod + theme(legend.position = "none"), 
          empty, key,
          labels = c('(a)', '(b)', '(c)', '(d)', '(e)', ' (f)', '', ''),
          label_size = 11,
          hjust = -4.5, vjust = 2.5,
          ncol = 2, align = "hv", axis = "l", rel_heights = c(1, 1, 1, 0.3)) -> ts_plot

## write plot
ggsave(filename = "./Figs/all_timeseries.png", device = "png", width = 11.5, height = 7, 
       dpi = 300) 


#Mean difference plot for Reviewer 1 response: 
ts %>% 
  filter(SIZESEX != "POP") %>%
  group_by(SIZESEX) %>%
  mutate(center_OCC = jtools::center(TEMP_OCC)) %>% #mean-center temp of occupancy 
  mutate(Diff = AVG_BT - TEMP_OCC) %>%
  mutate(center_Diff = AVG_BT - center_OCC) %>%
select(Year, SIZESEX, Diff, center_Diff)  %>%
  ggplot(aes(x = Year, y = Diff, group = SIZESEX, 
             shape = SIZESEX, color = SIZESEX))+
  geom_point()+
  geom_line()+
  scale_colour_manual(values = cbpalette, 
                      labels = c("Immature Females", "Mature Females", "Males 31 to 60 mm",
                                 "Males 61 to 90 mm", "Males 91 to 120 mm"))+
  scale_shape_manual(values = c(7:11),
                     labels = c("Immature Females", "Mature Females", "Males 31 to 60 mm",
                                "Males 61 to 90 mm", "Males 91 to 120 mm"))+
  labs(y = expression("Temperature Difference " ( degree~C)), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(axis.text.x  = element_text(size=11)) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 9.5),
        legend.title = element_blank())

## write plot
ggsave(filename = "./Figs/temp_diff.png", device = "png", width = 8, height = 6, 
       dpi = 300) 

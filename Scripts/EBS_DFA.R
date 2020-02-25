# notes ----
# Dynamic Factor Analysis to summarize environmental variation in cold pool COD/extent and bottom temps
# Mike Litzow
# last updated: 2020/2/5

# load ----
library(tidyverse)
library(MARSS)
library(broom)
library(cowplot)
library(FNGr) # for tickr function from Ben Williams, download package from GitHub

# set x axis labels
x_axis <- tickr(data = tibble(yr = 1988:2019), yr, 5)

# data ----
dat <- read.csv("./Data/Spatial_Env_TS.csv", row.names = 1)

# select the columns we want, and then transpose and turn into a matrix 
# (these are requirements for MARSS - column names need to be time steps)

dat <- dat %>%
  select(AVG_BT, CP_EXTENT, CP_COD)
head(dat)

# look at the distributions quickly!
plot.dat <- dat %>%
  gather()
  
ggplot(plot.dat, aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales="free") # not terrible?

# now set up for MARSS
dat <- as.matrix(t(dat))

# set up forms of R matrices - these are the four candidate error structures
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
# I'm using the default convergence criteria!

# make an object to save output
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  # allowing only one shared trend as we only have 3 time series!
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(dat, model=dfa.model,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data # unconstrained errors the best

# save the model comparison table
dir.create("output", showWarnings = FALSE)
write.csv(model.data, "./Output/environment dfa model selection.csv")

# now run the best model, save the trend estimates and CIs, and plot loadings/trend
model.list = list(A="zero", m=1, R="unconstrained")
mod <- MARSS(dat, model=model.list, z.score=TRUE, form="dfa")
summary(mod)


# get CI and plot loadings...
modCI <- MARSSparamCIs(mod)
plot.CI <- data.frame(names=rownames(dat), mean=modCI$par$Z, upCI=modCI$par.upCI$Z,
                      lowCI=modCI$par.lowCI$Z)
plot.CI$names <- reorder(plot.CI$names, plot.CI$mean)

# set colors and labels
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Envlabs <- c("Cold Pool Extent", "Bottom Temp", "Cold Pool COD")
dodge <- position_dodge(width=0.9)

#Plot Loadings 
plot.CI %>%
  mutate(names = c("Bottom Temp", "Cold Pool Extent", "Cold Pool COD")) %>%
  ggplot(aes(x = names, y = mean)) +
  geom_bar(position = "dodge", stat = "identity", fill = cb[2]) +
  geom_errorbar(aes(ymax = upCI, ymin = lowCI), position = dodge, width = 0.3) +
  geom_hline(yintercept = 0) +
  labs(y = "Loading", x = NULL) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) -> loadings
  

#Plot Shared Trend
trend.CI <- tidy(mod, type="states")
trend.CI$year <- 1988:2019

#Plot 
trend.CI %>%
  ggplot(aes(x = year, y = estimate)) +
  geom_line(color = cb[3], size = 1) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low), fill = cb[3], alpha = 0.2) +
  labs(y = "Shared Trend", x = NULL) +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  geom_hline(yintercept = 0) -> trend

#Combined plot for manuscript
plot_grid(trend, loadings, nrow = 1, rel_widths = c(.9, 0.5),
          labels = c('(a)', '(b)'),
          label_size = 11,
          hjust = -.2, vjust = 1.5)

## write plot
ggsave(filename = "./Figs/DFA.png", device = "png", width = 8, height =4, 
       dpi = 300) 

# Save the trend output
write.csv(trend.CI, "./Output/environment dfa trend.csv", row.names = F)
  #Trend Output was merged with "Spatial _Env_bySizeSex.csv" for LME/GLS modelling 

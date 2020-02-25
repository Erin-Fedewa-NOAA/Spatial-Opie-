# notes ----
#Erin Fedewa

#Goal: Mixed (LME) vrs fixed (GLS) effects models to test effects of environmental trend (via DFA analysis) on snow crab  
  #spatial extent (D95) and Lat center of dist (COD) after accounting for size/sex 
##See DFA script for DFA_Trend covariate development- trend is reduced from cold pool extent/center and avg temp 

#Approach (via Zuur- Mixed Effects)  
  #1) Center response variables to account for different intercepts by size/sex category
  #2) Fit full fixed effects model and compare correlation structures with AICc using REML
  #3) Re-fit lme random slope and gls full models with selected correlation structure and compare w/ AICc/REML
        #Does model selection favor a random effect for size/sex? 
  #4) D95 full fixed effects structure includes abundance as fixed effect
      # so re-fit with ML and compare nested models to evaluate support for fixed terms 
  #5) Refit selected model with REML


#load ----
library(tidyverse)
library(nlme)
library(MuMIn)
library(ggeffects)
library(cowplot)
library(forecast)
library(stargazer)

#data ----

dat <- read.csv("./Data/Spatial_Env_bySizeSex.csv")
head(dat)

dat %>%
  filter(SIZESEX != "POP") %>%
  group_by(SIZESEX) %>%
  mutate(center.D95 = jtools::center(D95)) %>%
  mutate(center.COD = jtools::center(LAT_COD)) ->dat #could also scale, but that would standardize variance 

# D95 Models ----

#1) Compare correlation structures for full Fixed Effects Model (Fixed: Trend and Abundance)
  corr0<-gls(center.D95~1+DFA_Trend+ABUN_mil_strat, data = dat) #OLS base model 
  acf(resid(corr0)) 
  corr1<-gls(center.D95~1+DFA_Trend+ABUN_mil_strat, correlation=corAR1(), data = dat) 

  #Specify Autocorrelation nested by Size/Sex (vrs default order of data by year) and compare
  corr2<-gls(center.D95~1+DFA_Trend+ABUN_mil_strat, correlation=corAR1(form = ~1|SIZESEX), data = dat) 
  AICc(corr0, corr1, corr2) #model 2 favored, but very marginal improvement
  
  #Check whether first order autocorrelation is sufficient 
  auto.arima(dat$center.D95, trace=1)
  #Looks like lag one is sufficient, lets go with corr2 structure 

#2)Re-fit lme random slope and gls full models with selected correlation structure and compare w/ AICc/REML
  #No need to look at random intercept  b/c response already scaled for different intercepts by size/sex
  
  #Full Fixed Effects Model --
  cm1<-gls(center.D95~1+DFA_Trend+ABUN_mil_strat, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat) 

  #Full Mixed Effects model with random slope on DFA trend for sizesex
  cm2<-lme(center.D95~1+DFA_Trend+ABUN_mil_strat, random= ~-1+DFA_Trend|SIZESEX, 
        method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
  summary(cm2)
  
  #Full Mixed Effects model with random slope on Abundance for sizesex
  cm3<-lme(center.D95~1+DFA_Trend+ABUN_mil_strat, random= ~-1+ABUN_mil_strat|SIZESEX, 
           method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
  summary(cm3)
  
  AICc(cm1, cm2, cm3) #Fixed effects model is supported 
    

#3) Test for best-supported fixed structure using likelihood ratio test comparing nested models 
 
  #Full model 
  cm5<-gls(center.D95~1+DFA_Trend+ABUN_mil_strat, method="ML", correlation=corAR1(form = ~1|SIZESEX), data = dat) 
  summary(cm5) #Abundance fixed effect not significant 

  #Drop abundance fixed effect
  cm6<-gls(center.D95~1+DFA_Trend, method="ML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
    AICc(cm5,cm6) #model without abundance is supported
  
#4) Re-fit final model with REML estimation 
  
    mod_final <-gls(center.D95~1+DFA_Trend, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
      summary(mod_final) #Env trend fixed effect significant
      
#Diagnostics
  plot(mod_final)
    resid<-resid(mod_final, type="normalized") #extract normalized residuals 
    fit<-fitted(mod_final)
  plot(fit, resid)
    abline(h=0)
     hist(resid) #Looks good
      acf(resid) #Also looks good 

#Extract results
  stargazer(mod_final, type = "text",
            digits = 3,
            star.cutoffs = c(0.05, 0.01, 0.001),
            digit.separator = "")

#Plots 
  #Color-blind palette
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  #Raw data plot with smoothed trend line 
dat %>%
  ggplot(aes(DFA_Trend, D95, col=SIZESEX, shape=SIZESEX)) +
  geom_point() +
  geom_smooth(method="lm", se=F) + #raw data smoothed
  labs(y="D95 (nm?)", x= "Environmental Trend") +
  scale_color_manual(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males"),
            values=cbPalette) +
   theme_bw() +theme(legend.title=element_blank())+
  theme(panel.grid = element_blank(),legend.position = "none")->g1


  # Extract the prediction data frame
  pred.mm <- ggpredict(mod_final, terms = c("DFA_Trend"))  # this gives overall population-level predictions for the model
  plot(pred.mm)
  
  #Fig 3 in MS: Plot the model predicted fit  
ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +    # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = dat,                      # adding the raw data (centered values)
               aes(x = DFA_Trend, y = center.D95, colour = SIZESEX)) + 
    labs(y=bquote('Mean-centered Snow Crab Areal Extent ('~nm^2~')'), x= "Shared Environmental Trend") +
    scale_color_manual(labels = c("Immature Females", "Mature Females", "Males 31 to 60mm", "Males 61 to 90mm","Males 91 to 120mm"),
                       values=cbPalette) +
    theme_bw() +theme(legend.title=element_blank()) +
    theme(panel.grid = element_blank())-> g2

ggsave(g2, filename = "./Figs/modelfit.png", device = "png", width = 7, height = 5, 
       dpi = 300)


# Center of Distribution Models ----

#1) Compare correlation structures for full Fixed Effects Model (Fixed: Trend)
corr0<-gls(center.COD~1+DFA_Trend, data = dat) #OLS base model 
  acf(resid(corr0)) 
corr1<-gls(center.COD~1+DFA_Trend, correlation=corAR1(), data = dat) 

#Specify Autocorrelation nested by Size/Sex (vrs default order of data by year) and compare
corr2<-gls(center.COD~1+DFA_Trend, correlation=corAR1(form = ~1|SIZESEX), data = dat) 
AICc(corr0, corr1, corr2) #model 1 AR form favored 

#Check whether first order autocorrelation is sufficient 
auto.arima(dat$center.COD, trace=1) #looks like ARMA(2,0,2) 
corr3<-gls(center.COD ~1+DFA_Trend, correlation = corARMA(p=2,q=2), data=dat)
AICc(corr1, corr3) #as indicated by auto.arima function, model 3 favored

#2)Re-fit lme random slope and gls full models with selected correlation structure and compare w/ AICc/REML

#Full Fixed Effects Model
mcentered1<-gls(center.COD~1+DFA_Trend, method="REML", correlation = corARMA(p=2,q=2), data = dat)

#Random slope for SIZESEX
mcentered2<-lme(center.COD~1+DFA_Trend, random= ~-1+DFA_Trend|SIZESEX, method="REML", correlation = corARMA(p=2,q=2), data = dat)
  AICc(mcentered1, mcentered2) #Model selection favors fixed effects model 

#3) Final model summary 
mod_final <- gls(center.COD~1+DFA_Trend, correlation = corARMA(p=2,q=2), data = dat)
  summary(mod_final)
  intervals(mod_final)
#DFA environmental trend fixed effect NOT significant

#Diagnostics
plot(mod_final)
  resid<-resid(mod_final, type="normalized") #extract normalized residuals 
  fit<-fitted(mod_final)
plot(fit, resid)
  hist(resid) 
  acf(resid) 

#Plots
  dat %>%
    ggplot(aes(DFA_Trend, LAT_COD, col=SIZESEX, shape=SIZESEX)) +
    geom_point() +
    geom_smooth(method="lm", se=F) + #raw data smoothed
    labs(y="Snow Crab Center of Distribution", x= "Environmental Trend") +
    scale_color_manual(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males"),
                       values=cbPalette) +
    theme_bw() +theme(legend.title=element_blank())+
    theme(panel.grid = element_blank(),legend.position = "none")->g3
  
# Extract the prediction data frame
  pred.mm <- ggpredict(mod_final, terms = c("DFA_Trend"))  # this gives overall predictions for the model
  pred.mm
  
# Plot the model predicted fit  
  ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +    # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = dat,                      # adding the raw data (centered values)
               aes(x = DFA_Trend, y = center.COD, colour = SIZESEX)) + 
    labs(y="Mean-centered Snow Crab \n Center of Distribution", x= "Shared Environmental Trend") +
    scale_color_manual(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males"),
                       values=cbPalette) +
    theme_bw() +theme(legend.title=element_blank()) +
    theme(panel.grid = element_blank()) -> g4
  

## combine plots
  plot_grid(g2 + theme(legend.position="none"), 
            g4 + theme(legend.position="none"),
            align = 'vh',
            labels = c("a)", "b)")) -> ts_plot
  
  # extract a legend that is laid out horizontally
  legend_b <- get_legend(g2 + 
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  
  # add the legend underneath the plot
  plot_grid(ts_plot, legend_b, ncol = 1, rel_heights = c(4, .2))
  
## write plot
  ggsave(filename = "./Output/modelfits.png", device = "png", width = 10, height = 5, 
         dpi = 300)
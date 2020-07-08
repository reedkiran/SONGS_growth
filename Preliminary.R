#SONGS PACL Data

##LOWSS (smooth fit- look for curve vs straight line)
##include standard length as covariate in model after potential transformation
##plot mean residuals graphically (plug into first graph)

#Packages
library(ggplot2)
library(car)
library(tidyverse)
library(plyr)

#Load Data
setwd("/Users/kiranreed/Desktop")
PACL_09_19 <- read.csv("PACL_data.csv")
SONGS_data <- read.csv("fish_collections_data_2009-2019.csv")
PACL_09_19_noNA<-drop_na(PACL_09_19,fish_growth_grams)

##########################################################################################

#Variance- Levene's Test: ASSUMPTION MET
PACL_09_19$year <- as.factor(PACL_09_19$year) #convert year to categorical
leveneTest(fish_growth_grams~reef_code, data=PACL_09_19)
boxplot(fish_growth_grams~reef_code, data=PACL_09_19) #visualize


#qqPlots
##no transformations
qqPlot(fish_growth_grams~reef_code, data=PACL_09_19)

##log transformed
PACL_09_19_noNA$log_growth<-log10(PACL_09_19_noNA$fish_growth_grams)
PACL_09_19_noNA$log_growth<-as.factor(PACL_09_19_noNA$log_growth)
PACL_09_19_noNA$log_growth<-revalue(PACL_09_19_noNA$log_growth, c("-Inf"="0.000000000000000000000001"))
PACL_09_19_noNA$log_growth<-as.numeric(PACL_09_19_noNA$log_growth)

qqPlot(log_growth~reef_code, data=PACL_09_19_noNA)

##square-root transformed
PACL_09_19_noNA$sqrt_growth<-sqrt(PACL_09_19_noNA$fish_growth_grams)

qqPlot(sqrt_growth~reef_code, data=PACL_09_19_noNA)

##########################################################################################

#ANOVA no covariates
aov.1 <- aov(fish_growth_grams~year*reef_code, data = PACL_09_19)
summary(aov.1)   

#ANOVA length covariate
aov.2 <- aov(fish_growth_grams~standard_length_mm+year*reef_code, data = PACL_09_19)
summary(aov.2)

##########################################################################################

#Standard Length by Growth Rate Plots
##LOWESS trend
ggplot(PACL_09_19, aes(x=standard_length_mm, y=fish_growth_grams)) +
  geom_point(alpha=0.5, aes(color=reef_code)) +
  labs(x="Standard Length (mm)", y= "Annual Fish Growth (grams)",
       title="PACL Growth Rate by Site (2009-2019)", legend="Site") +
  stat_smooth()+
  theme_classic()

##LOWESS trend for sqrt transformation
ggplot(PACL_09_19_noNA, aes(x=standard_length_mm, y=sqrt_growth)) +
  geom_point(alpha=0.5, aes(color=reef_code)) +
  labs(x="Standard Length (mm)", y= "Sqrt Annual Fish Growth (grams)",
       title="PACL Growth Rate by Site (2009-2019)", legend="Site") +
  stat_smooth()+
  theme_classic()

##Linear by Site trend line
ggplot(PACL_09_19, aes(x=standard_length_mm, y=fish_growth_grams, color=reef_code)) +
  geom_point(alpha=0.5) +
  labs(x="Standard Length (mm)", y= "Annual Fish Growth (grams)",
       title="PACL Growth Rate by Site (2009-2019)") +
  geom_smooth(method=lm,se=FALSE)+
  theme_classic()

##Linear trend Line
ggplot(PACL_09_19, aes(x=standard_length_mm, y=fish_growth_grams)) +
  geom_point(alpha=0.5, aes(color=reef_code)) +
  labs(x="Standard Length (mm)", y= "Annual Fish Growth (grams)",
       title="PACL Growth Rate by Site (2009-2019)", legend="Site") +
  geom_smooth(method=lm,se=FALSE)+
  theme_classic()

##Linear trend Line for sqrt transformed data
ggplot(PACL_09_19_noNA, aes(x=standard_length_mm, y=sqrt_growth)) +
  geom_point(alpha=0.5, aes(color=reef_code)) +
  labs(x="Standard Length (mm)", y= "Sqrt Annual Fish Growth (grams)",
       title="PACL Growth Rate by Site (2009-2019)", legend="Site") +
  geom_smooth(method=lm,se=FALSE)+
  theme_classic()

#Standardize growth by size with LM
growth.lm = lm(fish_growth_grams ~ standard_length_mm, data=PACL_09_19) 
summary(growth.lm)$r.squared
growth.res = resid(growth.lm)

PACL_09_19_noNA$resids<-growth.res

avg_site_yr.P<-PACL_09_19_noNA%>%
  group_by(reef_code,year)%>%
  summarize(avg=mean(fish_growth_grams,na.rm = FALSE),
            stdev=sd(fish_growth_grams,na.rm = FALSE),
            se=(sd(fish_growth_grams,na.rm = FALSE))/sqrt(NROW(fish_growth_grams)),
            resids=mean(resids))

avg_site_yr.P$year <- as.factor(avg_site_yr.P$year)

#Preliminary: average growth per year not standardized for size
ggplot(avg_site_yr.P, 
       aes(x=year, y=avg, group=reef_code, color=reef_code)) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se, width=0))+
  geom_line()+
  labs(x="Year", y= "Average Annual Fish Growth (grams)",
       title="PACL Growth Rate by Site (2009-2019)") +
  scale_colour_hue(name="Site")+
  theme_classic()

#Average growth per year standardized for size (residuals plot)
ggplot(avg_site_yr.P, 
       aes(x=year, y=resids, group=reef_code, color=reef_code)) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(ymin=resids-se, ymax=resids+se, width=0))+
  geom_line()+
  labs(x="Year", y= "Mean Residuals",
       title="PACL Growth Rate by Site (2009-2019)") +
  scale_colour_hue(name="Site")+
  theme_classic()

#Standardize growth by size with LM- sqrt transformed
growth.lm.T = lm(sqrt_growth ~ standard_length_mm, data=PACL_09_19_noNA) 
summary(growth.lm.T)$r.squared
growth.res.T = resid(growth.lm.T)

PACL_09_19_noNA$resids.T<-growth.res.T

avg_site_yr.P.T<-PACL_09_19_noNA%>%
  group_by(reef_code,year)%>%
  summarize(avg=mean(sqrt_growth,na.rm = FALSE),
            stdev=sd(sqrt_growth,na.rm = FALSE),
            se=(sd(sqrt_growth,na.rm = FALSE))/sqrt(NROW(fish_growth_grams)),
            resids=mean(resids.T))

avg_site_yr.P.T$year <- as.factor(avg_site_yr.P.T$year)

#Transformed average growth per year standardized for size (residuals plot)
ggplot(avg_site_yr.P.T, 
       aes(x=year, y=resids, group=reef_code, color=reef_code)) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(ymin=resids-se, ymax=resids+se, width=0))+
  geom_line()+
  labs(x="Year", y= "Mean Residuals",
       title="SQRT PACL Growth Rate by Site (2009-2019)") +
  scale_colour_hue(name="Site")+
  theme_classic()
##########################################################################################
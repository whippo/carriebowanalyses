###################################################################################
#                                                                                ##
# GLM models of comsumption assays at Carrie Bow Cay                             ##
# Data are current as of 2018-01-23                                              ##
# Data source: MarineGEO Program, Smithsonian                                    ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-01-23                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:
# This collection of GLMs explores the relationships between feeding pressure 
# (squidpops and weedpops) and fish community metrics at Carrie Bow Cay, Belize
# between 2015 and 2017.

# Required Files (check that script is loading latest version):
# 20171120_CBC_RLS_subsampled.csv
# Traits_all-species_edit.csv
# 20171122_squidpops_CBC.csv
# 20171120_weedpops.csv
# 20180124_CBC_biom.csv

# Associated Scripts:
# CBC_Manuscript_Figures.R
# RLS_M1M2_v1.R 

# TO DO

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# SUMMARY STATS                                                                   # #                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2018-01-23  Script Created

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(plyr)
library(tidyverse) # plots #data manipulation
library(ggpubr) # multi-plot visualizations
library(vegan) # diversity and MDS plots
library(splitstackshape) # data manipulation
library(viridis) # color palette
library(psych) # pairs analysis
library(lme4) # glm analyses
library(grid) # nMDS vectors
library(sjPlot) # quick summaries of models

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# set working directory
setwd("~/Dropbox (Personal)/Datasets/CBC_MS_2018")



############### RLS

# Import subsampled RLS survey data
RLS_subsample <- read.csv("20171120_CBC_RLS_subsampled.csv")
RLS_subsample$Year <- as.factor(RLS_subsample$Year)
RLS_subsample$Habitat <- factor(RLS_subsample$Habitat, levels = c("Fore Reef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))
glimpse(RLS_subsample)

# Import RLS biomass data

RLS_biom <- read.csv("20180124_CBC_biom.csv")



############### SQUIDPOPS

# read in csv (summarized data)
squidpops <- read.csv("20171122_squidpops_CBC.csv")

squidpops$Habitat <- factor(squidpops$Habitat, levels = c("fore reef", "patch reef", "mangrove", "seagrass", "sand"))

# change numbers to numeric
squidpops$Proportion.Missing.24.hours <- as.numeric(as.character(squidpops$Proportion.Missing.24.hours))
squidpops$Proportion.Missing.1.hour <- as.numeric(as.character(squidpops$Proportion.Missing.1.hour))

# change year to factor
squidpops$Year <- as.factor(squidpops$Year)


############### WEEDPOPS

# read in csv
weedpops <- read.csv("20171120_weedpops.csv")

# glimpse data
glimpse(weedpops)

weedpops$habitat <- factor(weedpops$habitat, levels = c("fore reef", "patch reef", "mangrove", "seagrass", "sand"))

# change year to factor
weedpops$year <- as.factor(weedpops$year)


############### JOIN ASSAY DATA WITH COMMUNITY DATA

# create dataframe of total biomass per site/time for RLS
biom_summary <- RLS_biom %>%
  group_by(Site.Name, Year, Habitat) %>%
  summarise(sum(Total.Biomass))
  # change column name
  names(biom_summary)[names(biom_summary)=="sum(Total.Biomass)"] <- "Total.Biomass"
  # log 10 the biomass
  biom_summary$Log.Biomass <- log10(biom_summary$Total.Biomass)
  # create site/time column to join with other data
  biom_summary <- biom_summary %>%
    unite(Site.Time, Site.Name, Year, sep="_", remove = FALSE)
  
# create site/time column for squidpop data joining
  squidpop_summ <- squidpops %>%
    select(Site, Year, Number.Deployed, Detachment.1.hour) %>%
    unite(Site.Time, Site, Year, sep ="_", remove = TRUE) 
  
# create site/time column for weedpop data joining and summarise
  weed_summ <- weedpops %>%
    group_by(site, year) %>%
    summarise(sum(number.recovered), sum(detachment.24hr)) 
  # change column name
  names(weed_summ)[names(weed_summ)=="sum(number.recovered)"] <- "Weed.Recovered"
  names(weed_summ)[names(weed_summ)=="sum(detachment.24hr)"] <- "Weed.Detached"
  # create site/time column to join with other data
  weed_summ <- weed_summ %>%
    unite(Site.Time, site, year, sep="_", remove = TRUE) %>%
    select(Site.Time, Weed.Recovered, Weed.Detached)
  
# create site/time column for RLS diversity data
  RLS_summ <- RLS_subsample
  RLS_summ$count <- 1
  RLS_summ <- RLS_summ %>%
    unite(Site.Time, Site.Name, Year, sep= "_", remove = TRUE) %>%
    group_by(Site.Time) %>%
    summarise(sum(count), sum(sum))
  #change column name
  names(RLS_summ)[names(RLS_summ)=="sum(count)"] <- "Raw.Richness"
  names(RLS_summ)[names(RLS_summ)=="sum(sum)"] <- "Fish.Abundance"
  RLS_summ$Log.Abundance <- log10(RLS_summ$Fish.Abundance)
  
GLM_data <- left_join(biom_summary, squidpop_summ, by = "Site.Time")  
GLM_data <- left_join(GLM_data, RLS_summ, by = "Site.Time")
GLM_data <- left_join(GLM_data, weed_summ, by = "Site.Time")

# make year a factor
GLM_data$Year <- as.factor(GLM_data$Year)


###################################################################################
# SUMMARY STATS                                                                   #
###################################################################################

############### SIMPLE PLOTS OF RELATIONSHIPS


# squidpops
plot(Detachment.1.hour/Number.Deployed ~ Raw.Richness, data = GLM_data)
plot(Detachment.1.hour/Number.Deployed ~ Log.Abundance, data = GLM_data)
plot(Detachment.1.hour/Number.Deployed ~ Log.Biomass, data = GLM_data)
 
# weedpops
plot(Weed.Detached/Weed.Recovered ~ Raw.Richness, data = GLM_data)
plot(Weed.Detached/Weed.Recovered ~ Log.Abundance, data = GLM_data)
plot(Weed.Detached/Weed.Recovered ~ Log.Biomass, data = GLM_data)


############### NORMALITY TESTS




############### GLM

# squidpops
SP_GLM <- glmer(formula = cbind(Detachment.1.hour, Number.Deployed - Detachment.1.hour) ~ Raw.Richness + Log.Abundance + Log.Biomass + Habitat + (1|Year) + (1|Site.Name), family = binomial(logit), data = GLM_data)

summary(SP_GLM)

# weedpops
WP_GLM <- glmer(formula = cbind(Weed.Detached, Weed.Recovered - Weed.Detached) ~ Raw.Richness + Log.Abundance + Log.Biomass + Habitat + (1|Year), family = binomial(logit), data = GLM_data)

summary(WP_GLM)

# change through time

BIOM_GLM <- lmer(formula = Log.Biomass ~ Habitat * Year + (1|Site.Name), data = GLM_data)
sjt.lmer(BIOM_GLM)

 ABUN_GLM <- lmer(formula = Log.Abundance ~ Habitat * Year + (1|Site.Name), data = GLM_data)
sjt.lmer(ABUN_GLM)

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#
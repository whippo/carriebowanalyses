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
# SUMMARY STATS                                                                   #
# LINEAR MODEL PLOTS                                                              #
# GLM                                                                             #
# BINOMIAL PLOT                                                                   #
# SEM                                                                             #
# 
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2018-02-06 SEM, GLM, and various plot labeled/added
# 2018-01-23  Script Created

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(tidyverse) # plots #data manipulation
library(ggpubr) # multi-plot visualizations
library(vegan) # diversity and MDS plots
library(splitstackshape) # data manipulation
library(viridis) # color palette
library(psych) # pairs analysis
library(lme4) # glm analyses
library(nlme) # SEM analysis
library(piecewiseSEM) # SEM
library(grid) # nMDS vectors
library(sjPlot) # quick summaries of models
library(piecewiseSEM) # computing pseudo-r for glms
library(DHARMa) # for residual diagnostics of glms (normality of residuals)
library(car) # test colinearity in GLM

############### vif.mer function

# find colinearity among variables

vif.mer <- function(fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v) ^ 0.5
  v <- diag(solve(v / (d %o% d)))
  names(v) <- nam
  v
}



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

# Import large fish dataset
RLS_large <- read.csv("20180124_CBC_largefish.csv")
RLS_large_summ <- RLS_large %>%
  group_by(Site.Name, Year) %>%
  summarise(sum(Total.Abundance))
names(RLS_large_summ)[names(RLS_large_summ) == "sum(Total.Abundance)"] <- "Large.Abun"
# create log abundance column
RLS_large_summ$Large.Log <- log10((RLS_large_summ$Large.Abun) + 1)

############### SQUIDPOPS

# read in csv (summarized data)
squidpops <- read.csv("20171122_squidpops_CBC.csv")

squidpops$Habitat <- factor(squidpops$Habitat, levels = c("fore reef", "patch reef", "mangrove", "seagrass", "sand"))

# change numbers to numeric
squidpops$Proportion.Missing.24.hours <- as.numeric(as.character(squidpops$Proportion.Missing.24.hours))
squidpops$Proportion.Missing.1.hour <- as.numeric(as.character(squidpops$Proportion.Missing.1.hour))

# change year to factor
squidpops$Year <- as.factor(squidpops$Year)

####### CODE ADAPTED FROM BRIAN CHENG TO EXPAND SQUIDPOP DETACHMENT TO
####### 0S AND 1S BINOMIAL.

# create column of retained squid
squidpops$Not.Detached <- squidpops$Number.Deployed - squidpops$Detachment.1.hour

# remove rows with NA

squidpops <- squidpops[!(is.na(squidpops$Detachment.1.hour)), ]

dataset <- squidpops # first clone your dataframe and call it 'dataset' so you don't have to rename stuff
# my columns were originally called "Alive" and "Dead", change these column names so you they are standardized at "success" (not detached) and "failure" (detached)
names(dataset)[names(dataset) == "Not.Detached"] <- c("success") # rename the dataframe columns
names(dataset)[names(dataset) == "Detachment.1.hour"] <- c("failure")
dataset.expanded <- dataset[0, ]
for (i in 1:length(dataset$success))
{
  if (dataset$success[i] > 0) {
    dataset.add.succ <- dataset[rep(i, dataset$success[i]), ]
    dataset.add.succ$success <- 1
    dataset.add.succ$failure <- 0
    dataset.expanded <- rbind(dataset.expanded, dataset.add.succ)
  }
  if (dataset$failure[i] > 0) {
    dataset.add.fail <- dataset[rep(i, dataset$failure[i]), ]
    dataset.add.fail$success <- 0
    dataset.add.fail$failure <- 1
    dataset.expanded <- rbind(dataset.expanded, dataset.add.fail)
  }
}

# Duplicate to clarify name
squidpops.expanded <- dataset.expanded

############### WEEDPOPS

# read in csv
weedpops <- read.csv("20171120_weedpops.csv")

# glimpse data
glimpse(weedpops)

weedpops$habitat <- factor(weedpops$habitat, levels = c("fore reef", "patch reef", "mangrove", "seagrass", "sand"))

# change year to factor
weedpops$year <- as.factor(weedpops$year)

####### CODE ADAPTED FROM BRIAN CHENG TO EXPAND WEEDPOP DETACHMENT TO
####### 0S AND 1S BINOMIAL.

# create column of retained squid
weedpops$Not.Detached <- weedpops$number.recovered - weedpops$detachment.24hr

# remove rows with NA

weedpops <- weedpops[!(is.na(weedpops$detachment.24hr)), ]

dataset <- weedpops # first clone your dataframe and call it 'dataset' so you don't have to rename stuff
# my columns were originally called "Alive" and "Dead", change these column names so you they are standardized at "success" (not detached) and "failure" (detached)
names(dataset)[names(dataset) == "Not.Detached"] <- c("success") # rename the dataframe columns
names(dataset)[names(dataset) == "detachment.24hr"] <- c("failure")
dataset.expanded <- dataset[0, ]
for (i in 1:length(dataset$success))
{
  if (dataset$success[i] > 0) {
    dataset.add.succ <- dataset[rep(i, dataset$success[i]), ]
    dataset.add.succ$success <- 1
    dataset.add.succ$failure <- 0
    dataset.expanded <- rbind(dataset.expanded, dataset.add.succ)
  }
  if (dataset$failure[i] > 0) {
    dataset.add.fail <- dataset[rep(i, dataset$failure[i]), ]
    dataset.add.fail$success <- 0
    dataset.add.fail$failure <- 1
    dataset.expanded <- rbind(dataset.expanded, dataset.add.fail)
  }
}

# Duplicate to clarify name
weedpops.expanded <- dataset.expanded


############### JOIN ASSAY DATA WITH COMMUNITY DATA

# create dataframe of total biomass per site/time for RLS
biom_summary <- RLS_biom %>%
  group_by(Site.Name, Year, Habitat) %>%
  summarise(sum(Total.Biomass))
# change column name
names(biom_summary)[names(biom_summary) == "sum(Total.Biomass)"] <- "Total.Biomass"
# log 10 the biomass
biom_summary$Log.Biomass <- log10(biom_summary$Total.Biomass)
# create site/time column to join with other data
biom_summary <- biom_summary %>%
  unite(Site.Time, Site.Name, Year, sep = "_", remove = FALSE)

# create site/time for large fish abundance
RLS_large_summ <- RLS_large_summ %>%
  select(Site.Name, Year, Large.Log) %>%
  unite(Site.Time, Site.Name, Year, sep = "_", remove = TRUE)

# create site/time column for squidpop data joining
squidpop_summ <- squidpops.expanded %>%
  select(Site, Year, failure) %>%
  unite(Site.Time, Site, Year, sep = "_", remove = TRUE)
# change name of failure to Squidpop.Detached (0=no 1=yes)
names(squidpop_summ)[names(squidpop_summ) == "failure"] <- "Squidpop.Detached"

# create site/time column for weedpop data joining and summarise
weedpop_summ <- weedpops.expanded %>%
  select(site, year, failure) %>%
  unite(Site.Time, site, year, sep = "_", remove = TRUE)
names(weedpop_summ)[names(weedpop_summ) == "failure"] <- "Weedpop.Detached"

# create site/time column for RLS diversity data
RLS_summ <- RLS_subsample
RLS_summ$count <- 1
RLS_summ <- RLS_summ %>%
  unite(Site.Time, Site.Name, Year, sep = "_", remove = TRUE) %>%
  group_by(Site.Time) %>%
  summarise(sum(count), sum(sum))
# change column name
names(RLS_summ)[names(RLS_summ) == "sum(count)"] <- "Raw.Richness"
names(RLS_summ)[names(RLS_summ) == "sum(sum)"] <- "Fish.Abundance"
RLS_summ$Log.Abundance <- log10(RLS_summ$Fish.Abundance)

# Join all for squidpop analysis
GLM_data_SP <- left_join(biom_summary, squidpop_summ, by = "Site.Time")
GLM_data_SP <- left_join(GLM_data_SP, RLS_summ, by = "Site.Time")
GLM_data_SP <- left_join(GLM_data_SP, RLS_large_summ, by = "Site.Time")

# make year a factor
GLM_data_SP$Year <- as.factor(GLM_data_SP$Year)

# Join all for weedpop analysis
GLM_data_WP <- left_join(biom_summary, weedpop_summ, by = "Site.Time")
GLM_data_WP <- left_join(GLM_data_WP, RLS_summ, by = "Site.Time")
GLM_data_WP <- left_join(GLM_data_WP, RLS_large_summ, by = "Site.Time")

# make year a factor
GLM_data_WP$Year <- as.factor(GLM_data_WP$Year) 

#write.csv(GLM_data_SP, "GLM_data_SP.csv")





###################################################################################
# LINEAR MODEL PLOTS                                                              #
###################################################################################

############### SIMPLE PLOTS OF RELATIONSHIPS


# squidpops
plot(Squidpop.Detached ~ Raw.Richness, data = GLM_data_SP)
plot(Squidpop.Detached ~ Log.Abundance, data = GLM_data_SP)
plot(Squidpop.Detached ~ Log.Biomass, data = GLM_data_SP)

# weedpops NEED TO UPDATE TO BINOMIAL
plot(Weed.Detached / Weed.Recovered ~ Raw.Richness, data = GLM_data)
plot(Weed.Detached / Weed.Recovered ~ Log.Abundance, data = GLM_data)
plot(Weed.Detached / Weed.Recovered ~ Log.Biomass, data = GLM_data)





###################################################################################
# GLM                                                                             #
###################################################################################

# Notes with Jon.

# cbind for column successs/colum failure. Should parse out to 0 and 1. Expand grid? Brian Cheng has code?

# Abundance and biomass are highly correlated. Model can't distribute because of this phenomena. Creates huge standard error, can't tell diff between these things.

# Jon's worked focused on abundance but did bite rate per biomass?

# car package has vif function, computes degree of colinearity between variables. Allows you to choose one varible or the other if they are correlate.

# richness might be correlated too.

# (1|Site.Name/Year) use this syntax, acknowleges resampled sites through time. If you get large SD associated with site, model accounted for lots of varition within site. If number is small, very little variation, may want to remove? Or could do psuedo-r square. piecewiseSEM function rsquares -> marginal r squared, output similar to lm and other measure that accounts for random factors. conditional r tell how much is accounted for my random effects.

# colinearity is major issue, start of model selection. Model selection can be self defeating, can't talk about the things that are interesting and mechanistic. Not a fishing expedition.

# normality test. residuals have to be normal, not predictors. Settle on model first, then check normality. Use Dharma package to check normality. Center and divide by SD of variables to scale.

# You know what's going on out there, you know the system. See if the model makes sense. Intuition should trump everything.

# expand.grid <- create 0 and 1 from proportions?



################ HABITAT AS A PREDICTOR

# squidpops glm, binomial model tested against fish metrics and habitat type, with
# year and site as random variables.
SP_GLM <- glmer(formula = Squidpop.Detached ~ Raw.Richness + Log.Biomass + Log.Abundance + Large.Log + Habitat + (1 | Site.Name / Year), family = binomial(logit), data = GLM_data_SP)
sjt.glmer(SP_GLM)
summary(SP_GLM)

# obtaining pseudo-r values for the model
rsquared(SP_GLM)

# run DHARMa simulation for residuals (testing normality)
simulationOutput <- simulateResiduals(fittedModel = SP_GLM, n = 250)

# visualize the output

plotSimulatedResiduals(simulationOutput = simulationOutput)

# run colinearity analysis for factors in analaysis (vif.mer function above)
vif.mer(SP_GLM)




# weedpops glm, binomial model tested against fish metrics and habitat type, with
# year and site as random variables.
WP_GLM <- glmer(formula = Weedpop.Detached ~ Raw.Richness + Log.Biomass + Log.Abundance + Large.Log + Habitat + (1 | Site.Name / Year), family = binomial(logit), data = GLM_data_WP)
sjt.glmer(WP_GLM)
summary(WP_GLM)

# obtaining pseudo-r values for the model
rsquared(WP_GLM)

# run DHARMa simulation for residuals (testing normality)
simulationOutput <- simulateResiduals(fittedModel = WP_GLM, n = 250)

# visualize the output

plotSimulatedResiduals(simulationOutput = simulationOutput)

# run colinearity analysis for factors in analaysis (vif.mer function above)
vif.mer(WP_GLM)



################ HABITAT AS A RANDOM VARIABLE


SP_GLM_habrand <- glmer(formula = Squidpop.Detached ~ Raw.Richness + Log.Biomass + Large.Log + (1 | Site.Name / Habitat / Year), family = binomial(logit), data = GLM_data_SP)
sjt.glmer(SP_GLM_habrand)
summary(SP_GLM_habrand)

# obtaining pseudo-r values for the model
rsquared(SP_GLM_habrand)

# run DHARMa simulation for residuals (testing normality)
simulationOutput <- simulateResiduals(fittedModel = SP_GLM_habrand, n = 250)

# visualize the output

plotSimulatedResiduals(simulationOutput = simulationOutput)

# run colinearity analysis for factors in analaysis (vif.mer function above)
vif.mer(SP_GLM_habrand)



###################################################################################
# BINOMIAL PLOT                                                                   #
###################################################################################


# VISUALIZATION CURVE

SP_GLM_habrand_curve <- glm(formula = Squidpop.Detached ~ Raw.Richness, family = binomial(logit), data = GLM_data_SP)

plot(Squidpop.Detached ~ Raw.Richness, data = GLM_data_SP)

curve(predict(SP_GLM_habrand_curve,data.frame(Raw.Richness=x),type="resp"),add=TRUE) 

# obtaining pseudo-r values for the model
rsquared(SP_GLM_habrand)

# run DHARMa simulation for residuals (testing normality)
simulationOutput <- simulateResiduals(fittedModel = SP_GLM_habrand, n = 250)

# visualize the output

plotSimulatedResiduals(simulationOutput = simulationOutput)

# run colinearity analysis for factors in analaysis (vif.mer function above)
vif.mer(SP_GLM_habrand)










# change through time

BIOM_GLM <- lmer(formula = Log.Biomass ~ Habitat * Year + (1 | Site.Name), data = GLM_data)
sjt.lmer(BIOM_GLM)

ABUN_GLM <- lmer(formula = Log.Abundance ~ Habitat * Year + (1 | Site.Name), data = GLM_data)
sjt.lmer(ABUN_GLM)



###################################################################################
# SEM                                                                             #
###################################################################################
# from Jon Lefcheck

### SEM for MarineGEO

# Load required libraries
#library(lme4)
#library(nlme)
#devtools::install_github("jslefche/piecewiseSEM@2.0")
#library(piecewiseSEM) # welcome message indicating version

# Load data
data <- read.csv("GLM_data_SP.csv")

data$Log.Richness <- log10(data$Raw.Richness)

# Create models going into SEM
model <- psem(
  
  diversity = lme(Log.Richness ~ Habitat ,random = ~ 1 | Year, data),
  
  biomass = lme(Total.Biomass ~ Habitat, random = ~ 1 | Year,  data),
  
  predation = glmer(Squidpop.Detached ~ Habitat + Total.Biomass + Log.Richness + (1 | Site.Name / Year), family = binomial, data), # logit link is the default
  
  Log.Richness %~~% Total.Biomass, # correlated error
  
  data = data
)

summary(model)




#####
# <<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#



####### SCRATCH PAD










## 1. center and scale predictors:
ss.SP <- transform(GLM_data_SP, Log.Biomass = scale(Log.Biomass), Raw.Richness = scale(Raw.Richness), Log.Abundance = scale(Log.Abundance), Large.Log = scale(Large.Log))
SP_GLM_scaled <- update(SP_GLM, data = ss.SP)


diag.vals <- getME(WP_GLM, "theta")[getME(WP_GLM, "lower") == 0]
any(diag.vals < 1e-6) # FALSE
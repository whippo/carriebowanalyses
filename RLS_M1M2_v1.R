###################################################################################
#                                                                                ##
# Reef Life Survey: Analyses of fish communities (all methods)                   ##
# Data are current as of 20171120                                                ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 20171120                                                          ##
#                                                                                ##
###################################################################################

# The aim of this script is to standardize the M1 and M2 data from each RLS survey
# so that a global richness and diversity value can be found for each site/time
# combination, without having to split it among survey types. 

# TO DO

# 2017 was submitted with different site names. Script was adjusted to deal with 
# these differences. Need to harmonize names across all years and make 2017 script
# same as previous years.

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# COMBINED M1/M2                                                                  #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 20171119 Corrected mislabeled Seagrass Habitat (was listed 'Patch Reef')
# 20171113 Added 2017 data scripts
# 20170413 File created from RLS_Script_v5.R, separated out M1/M2 combination analysis

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(nlme)
library(plyr)
library(car)
library(ggplot2)
library(reshape2) 
library(dplyr)
library(tidyr)
library(rfishbase)
library(vegan)
library(grid)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# Before import, be sure to change first column header "ID" to "Line_ID" in .csv. 
# Delete lat/long columns, and first empty row if present. Add Habitat Column.

# Change working directory
setwd("~/Dropbox (Personal)/R_Scripts/CBC 2017/RLS")

# Import the survey data (version 2 has correct Block values)
RLS_survey_data <- read.csv("CBCRLS_CORE.csv")

glimpse(RLS_survey_data)

# Create unique ID for each survey
RLS_survey_data <- RLS_survey_data %>%
  unite_("event", c("Site.Name", "Date", "Block"), sep = "_", remove = FALSE)

# make unique ID a factor
RLS_survey_data$event <- as.factor(RLS_survey_data$event)

# replace NAs with 0 
RLS_survey_data[is.na(RLS_survey_data)] <- 0

#remove method 1 with no given abundance
RLS_survey_data <- RLS_survey_data %>%
  subset(Total !="0")

# fix incorrect taxa
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Sphoeroides spenglerii"] <- "Sphoeroides spengleri"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Xyricthys splendens"] <- "Xyrichtys splendens"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Ctenogonius stigmaturus"] <- "Ctenogobius stigmaturus"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Ophistignathus sp."] <- "Opistognathus sp."
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Ctenogonius stigmaturus"] <- "Ctenogobius stigmaturus"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Hyppolytidae"] <- "Hippolytidae"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Paradiplogrammus bairdi"] <- "Callionymus bairdi"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Cryptotoums roseus"] <- "Cryptotomus roseus"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Eucinostoums gula"] <- "Eucinostomus gula"
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Atherinidae"] <- "Atherinid spp."
levels(RLS_survey_data$Species)[levels(RLS_survey_data$Species)=="Clupeoid spp."] <- "Clupeidae"


# remove unknown and freshwater species
RLS_survey_data <- RLS_survey_data %>%
  filter(Species != c("Gambusia sp.")) %>%
  filter(Species != c("Scarine spp."))

# all fish community all years 
RLS_1_2 <- RLS_survey_data %>%
  filter(Inverts < 1)
# remove incidental inverts
RLS_1_2 <- RLS_1_2 %>%
  filter(Species != "Hippolytidae") %>%
  filter(Species != "Mithracid spp.")%>%
  filter(Species != "Strombus gigas") %>%
  filter(Species != "Sepioteuthis sepioidea")




###################################################################################
# COMBINED M1/M2                                                                  #
###################################################################################

M1M2 <- RLS_1_2

# SUBSETTING M1 AND COMBINING WITH M2 FOR NEW SUBSAMPLED GLOBAL DATASET

# Calculate richness for all M2 surveys
M2 <- M1M2 %>%
  filter(Method == "2")

detach(package:plyr)
M2rich <- M2 %>%
  group_by(Site.Name, Habitat, Year, Species) %>%
  summarise(sum(Total))

# create M2 species list for each site/year so they can be targeted for removal in M1
BGma15 <- M2rich %>%
  filter(Site.Name == "Blueground Mangrove") %>%
  filter(Year == "2015")
CBpa15 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Patch Reef") %>%
  filter(Year == "2015")
CBfr15 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Reef") %>%
  filter(Year == "2015")
CBsa15 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Sand") %>%
  filter(Year == "2015")
CBsg15 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Seagrass") %>%
  filter(Year == "2015")
CUpr15 <- M2rich %>%
  filter(Site.Name == "Curlew Patch Reef") %>%
  filter(Year == "2015")
CUsa15 <- M2rich %>%
  filter(Site.Name == "Curlew Sand") %>%
  filter(Year == "2015")
CUsg15 <- M2rich %>%
  filter(Site.Name == "Curlew Seagrass") %>%
  filter(Year == "2015")
HRpr15 <- M2rich %>%
  filter(Site.Name == "House Reef") %>%
  filter(Year == "2015")
SRfr15 <- M2rich %>%
  filter(Site.Name == "South Reef") %>%
  filter(Year == "2015")
TBma15 <- M2rich %>%
  filter(Site.Name == "Tobacco Mangrove") %>%
  filter(Year == "2015")
TBfr15 <- M2rich %>%
  filter(Site.Name == "Tobacco Reef") %>%
  filter(Year == "2015")
TBsa15 <- M2rich %>%
  filter(Site.Name == "Tobacco Sand") %>%
  filter(Year == "2015")



CBpa16 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Patch Reef") %>%
  filter(Year == "2016")
CBfr16 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Reef") %>%
  filter(Year == "2016")
CBsa16 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Sand") %>%
  filter(Year == "2016")
CBsg16 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Seagrass") %>%
  filter(Year == "2016")
CUpr16 <- M2rich %>%
  filter(Site.Name == "Curlew Patch Reef") %>%
  filter(Year == "2016")
CUsa16 <- M2rich %>%
  filter(Site.Name == "Curlew Sand") %>%
  filter(Year == "2016")
CUsg16 <- M2rich %>%
  filter(Site.Name == "Curlew Seagrass") %>%
  filter(Year == "2016")
HRpr16 <- M2rich %>%
  filter(Site.Name == "House Reef") %>%
  filter(Year == "2016")
SRfr16 <- M2rich %>%
  filter(Site.Name == "South Reef") %>%
  filter(Year == "2016")
TBma16 <- M2rich %>%
  filter(Site.Name == "Tobacco Mangrove") %>%
  filter(Year == "2016")
TBfr16 <- M2rich %>%
  filter(Site.Name == "Tobacco Reef") %>%
  filter(Year == "2016")
TBsa16 <- M2rich %>%
  filter(Site.Name == "Tobacco Sand") %>%
  filter(Year == "2016")


BGma17 <- M2rich %>%
  filter(Site.Name == "Blueground Mangrove") %>%
  filter(Year == "2017")
CBpa17 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Patch Reef") %>%
  filter(Year == "2017")
CBfr17 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Reef") %>%
  filter(Year == "2017")
CBsa17 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Sand") %>%
  filter(Year == "2017")
CBsg17 <- M2rich %>%
  filter(Site.Name == "Carrie Bow Seagrass") %>%
  filter(Year == "2017")
CUpr17 <- M2rich %>%
  filter(Site.Name == "Curlew Patch Reef") %>%
  filter(Year == "2017")
CUsa17 <- M2rich %>%
  filter(Site.Name == "Curlew Sand") %>%
  filter(Year == "2017")
CUsg17 <- M2rich %>%
  filter(Site.Name == "Curlew Seagrass") %>%
  filter(Year == "2017")
HRpr17 <- M2rich %>%
  filter(Site.Name == "House Reef") %>%
  filter(Year == "2017")
SRfr17 <- M2rich %>%
  filter(Site.Name == "South Reef") %>%
  filter(Year == "2017")
TBma17 <- M2rich %>%
  filter(Site.Name == "Tobacco Mangrove") %>%
  filter(Year == "2017")
TBfr17 <- M2rich %>%
  filter(Site.Name == "Tobacco Reef") %>%
  filter(Year == "2017")
TBsa17 <- M2rich %>%
  filter(Site.Name == "Tobacco Sand") %>%
  filter(Year == "2017")
TCma17 <- M2rich %>%
  filter(Site.Name == "Twin Cays Mangrove") %>%
  filter(Year == "2017") 
TCsg17 <- M2rich %>%
  filter(Site.Name == "Twin Cays Seagrass") %>%
  filter(Year == "2017")
TBsg17 <- M2rich %>%
  filter(Site.Name == "Tobacco Seagrass") %>%
  filter(Year == "2017")



# Create M1/0 df per site/year
M10 <- M1M2 %>%
  filter(Method != "2")

detach(package:plyr)
M10rich <- M10 %>%
  group_by(Site.Name, Habitat, Year, Species) %>%
  summarise(sum(Total))

# create M10 species list for each site/year and remove redundant M2s
BGma15_10 <- M10rich %>%
  filter(Site.Name == "Blueground Mangrove") %>%
  filter(Year == "2015")
CBpa15_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Patch Reef") %>%
  filter(Year == "2015")
CBfr15_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Reef") %>%
  filter(Year == "2015")
CBsa15_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Sand") %>%
  filter(Year == "2015")
CBsg15_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Seagrass") %>%
  filter(Year == "2015")
CUpr15_10 <- M10rich %>%
  filter(Site.Name == "Curlew Patch Reef") %>%
  filter(Year == "2015")
CUsa15_10 <- M10rich %>%
  filter(Site.Name == "Curlew Sand") %>%
  filter(Year == "2015")
CUsg15_10 <- M10rich %>%
  filter(Site.Name == "Curlew Seagrass") %>%
  filter(Year == "2015")
HRpr15_10 <- M10rich %>%
  filter(Site.Name == "House Reef") %>%
  filter(Year == "2015")
SRfr15_10 <- M10rich %>%
  filter(Site.Name == "South Reef") %>%
  filter(Year == "2015")
TBma15_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Mangrove") %>%
  filter(Year == "2015")
TBfr15_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Reef") %>%
  filter(Year == "2015")
TBsa15_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Sand") %>%
  filter(Year == "2015")
TCma15_10 <- M10rich %>%
  filter(Site.Name == "Twin Cays Mangrove") %>%
  filter(Year == "2015")
TCsg15_10 <- M10rich %>%
  filter(Site.Name == "Twin Cays Seagrass") %>%
  filter(Year == "2015")


BGma16_10 <- M10rich %>%
  filter(Site.Name == "Blueground Mangrove") %>%
  filter(Year == "2016")
CBpa16_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Patch Reef") %>%
  filter(Year == "2016")
CBfr16_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Reef") %>%
  filter(Year == "2016")
CBsa16_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Sand") %>%
  filter(Year == "2016")
CBsg16_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Seagrass") %>%
  filter(Year == "2016")
CUpr16_10 <- M10rich %>%
  filter(Site.Name == "Curlew Patch Reef") %>%
  filter(Year == "2016")
CUsa16_10 <- M10rich %>%
  filter(Site.Name == "Curlew Sand") %>%
  filter(Year == "2016")
CUsg16_10 <- M10rich %>%
  filter(Site.Name == "Curlew Seagrass") %>%
  filter(Year == "2016")
HRpr16_10 <- M10rich %>%
  filter(Site.Name == "House Reef") %>%
  filter(Year == "2016")
SRfr16_10 <- M10rich %>%
  filter(Site.Name == "South Reef") %>%
  filter(Year == "2016")
TBma16_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Mangrove") %>%
  filter(Year == "2016")
TBfr16_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Reef") %>%
  filter(Year == "2016")
TBsa16_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Sand") %>%
  filter(Year == "2016")
TCma16_10 <- M10rich %>%
  filter(Site.Name == "Twin Cays Mangrove") %>%
  filter(Year == "2016")
TCsg16_10 <- M10rich %>%
  filter(Site.Name == "Twin Cays Seagrass") %>%
  filter(Year == "2016")

BGma17_10 <- M10rich %>%
  filter(Site.Name == "Blueground Mangrove") %>%
  filter(Year == "2017")
CBpa17_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Patch Reef") %>%
  filter(Year == "2017")
CBfr17_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Reef") %>%
  filter(Year == "2017")
CBsa17_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Sand") %>%
  filter(Year == "2017")
CBsg17_10 <- M10rich %>%
  filter(Site.Name == "Carrie Bow Seagrass") %>%
  filter(Year == "2017")
CUpr17_10 <- M10rich %>%
  filter(Site.Name == "Curlew Patch Reef") %>%
  filter(Year == "2017")
CUsa17_10 <- M10rich %>%
  filter(Site.Name == "Curlew Sand") %>%
  filter(Year == "2017")
CUsg17_10 <- M10rich %>%
  filter(Site.Name == "Curlew Seagrass") %>%
  filter(Year == "2017")
HRpr17_10 <- M10rich %>%
  filter(Site.Name ==  "House Reef") %>%
  filter(Year == "2017")
SRfr17_10 <- M10rich %>%
  filter(Site.Name == "South Reef") %>%
  filter(Year == "2017")
TBma17_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Mangrove") %>%
  filter(Year == "2017")
TBfr17_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Reef") %>%
  filter(Year == "2017")
TBsa17_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Sand") %>%
  filter(Year == "2017")
TCma17_10 <- M10rich %>%
  filter(Site.Name == "Twin Cays Mangrove") %>%
  filter(Year == "2017")
TCsg17_10 <- M10rich %>%
  filter(Site.Name == "Twin Cays Seagrass") %>%
  filter(Year == "2017")
TBsg17_10 <- M10rich %>%
  filter(Site.Name == "Tobacco Seagrass") %>%
  filter(Year == "2017")

# remove matching M2s from M1 df and subsample 1 in 5

########
# 2015 #
########

# seeds <- sample(1:1000, 50, replace = FALSE)
detach(package:plyr)

#BGma15 no matches
# create pool of M1 only fishes to draw from
names(BGma15_10)[names(BGma15_10)=="sum(Total)"] <- "sum"
BGma15_10_pool <- rep(BGma15_10$Species, BGma15_10$sum)
BGma15_10_pool <- factor(BGma15_10_pool)
set.seed(332)
BGma15_sampled <- as.data.frame(table((sample(BGma15_10_pool, sum(BGma15_10$sum)/5, replace = FALSE, prob = rep(1/sum(BGma15_10$sum), length(BGma15_10_pool))))))
BGma15_sampled$Site.Name <- "Blueground Mangrove"
BGma15_sampled$Habitat <- "Mangrove"  
BGma15_sampled$Year <- "2015"
names(BGma15_sampled)[names(BGma15_sampled)=="Var1"] <- "Species"
names(BGma15_sampled)[names(BGma15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(BGma15)[names(BGma15)=="sum(Total)"] <- "sum"
BGma15$Year <- as.character(BGma15$Year)
BGma15_final <- bind_rows(BGma15, BGma15_sampled)
BGma15_final <- BGma15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBpa15_remove <- CBpa15$Species[(CBpa15$Species %in% CBpa15_10$Species)]
CBpa15_remove <- factor(CBpa15_remove)
CBpa15_10 <- CBpa15_10 %>%
  filter(!(Species %in% CBpa15_remove))
names(CBpa15_10)[names(CBpa15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBpa15_10_pool <- rep(CBpa15_10$Species, CBpa15_10$sum)
CBpa15_10_pool <- factor(CBpa15_10_pool)
set.seed(945)
CBpa15_sampled <- as.data.frame(table((sample(CBpa15_10_pool, sum(CBpa15_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBpa15_10$sum), length(CBpa15_10_pool))))))
CBpa15_sampled$Site.Name <- "Carrie Bow Patch Reef"
CBpa15_sampled$Habitat <- "Patch Reef"  
CBpa15_sampled$Year <- "2015"
names(CBpa15_sampled)[names(CBpa15_sampled)=="Var1"] <- "Species"
names(CBpa15_sampled)[names(CBpa15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBpa15)[names(CBpa15)=="sum(Total)"] <- "sum"
CBpa15$Year <- as.character(CBpa15$Year)
CBpa15_final <- bind_rows(CBpa15, CBpa15_sampled)
CBpa15_final <- CBpa15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBfr15_remove <- CBfr15$Species[(CBfr15$Species %in% CBfr15_10$Species)]
CBfr15_remove <- factor(CBfr15_remove)
# check for actual overlap:
CBfr15_remove
CBfr15_10 <- CBfr15_10 %>%
  filter(!(Species %in% CBfr15_remove))
names(CBfr15_10)[names(CBfr15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBfr15_10_pool <- rep(CBfr15_10$Species, CBfr15_10$sum)
CBfr15_10_pool <- factor(CBfr15_10_pool)
set.seed(178)
CBfr15_sampled <- as.data.frame(table((sample(CBfr15_10_pool, sum(CBfr15_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBfr15_10$sum), length(CBfr15_10_pool))))))
CBfr15_sampled$Site.Name <- "Carrie Bow Reef"
CBfr15_sampled$Habitat <- "Fore Reef"  
CBfr15_sampled$Year <- "2015"
names(CBfr15_sampled)[names(CBfr15_sampled)=="Var1"] <- "Species"
names(CBfr15_sampled)[names(CBfr15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBfr15)[names(CBfr15)=="sum(Total)"] <- "sum"
CBfr15$Year <- as.character(CBfr15$Year)
CBfr15_final <- bind_rows(CBfr15, CBfr15_sampled)
CBfr15_final <- CBfr15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBsa15_remove <- CBsa15$Species[(CBsa15$Species %in% CBsa15_10$Species)]
CBsa15_remove <- factor(CBsa15_remove)
# check for actual overlap:
CBsa15_remove
CBsa15_10 <- CBsa15_10 %>%
  filter(!(Species %in% CBsa15_remove))
names(CBsa15_10)[names(CBsa15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw saom
CBsa15_10_pool <- rep(CBsa15_10$Species, CBsa15_10$sum)
CBsa15_10_pool <- factor(CBsa15_10_pool)
set.seed(930)
CBsa15_sampled <- as.data.frame(table((sample(CBsa15_10_pool, sum(CBsa15_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBsa15_10$sum), length(CBsa15_10_pool))))))
CBsa15_sampled$Site.Name <- "Carrie Bow Sand"
CBsa15_sampled$Habitat <- "Sand"  
CBsa15_sampled$Year <- "2015"
names(CBsa15_sampled)[names(CBsa15_sampled)=="Var1"] <- "Species"
names(CBsa15_sampled)[names(CBsa15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBsa15)[names(CBsa15)=="sum(Total)"] <- "sum"
CBsa15$Year <- as.character(CBsa15$Year)
CBsa15_final <- bind_rows(CBsa15, CBsa15_sampled)
CBsa15_final <- CBsa15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CUpr15_remove <- CUpr15$Species[(CUpr15$Species %in% CUpr15_10$Species)]
CUpr15_remove <- factor(CUpr15_remove)
# check for actual overlap:
CUpr15_remove
CUpr15_10 <- CUpr15_10 %>%
  filter(!(Species %in% CUpr15_remove))
names(CUpr15_10)[names(CUpr15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUpr15_10_pool <- rep(CUpr15_10$Species, CUpr15_10$sum)
CUpr15_10_pool <- factor(CUpr15_10_pool)
set.seed(569)
CUpr15_sampled <- as.data.frame(table((sample(CUpr15_10_pool, sum(CUpr15_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUpr15_10$sum), length(CUpr15_10_pool))))))
CUpr15_sampled$Site.Name <- "Curlew Patch Reef"
CUpr15_sampled$Habitat <- "Patch Reef"  
CUpr15_sampled$Year <- "2015"
names(CUpr15_sampled)[names(CUpr15_sampled)=="Var1"] <- "Species"
names(CUpr15_sampled)[names(CUpr15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUpr15)[names(CUpr15)=="sum(Total)"] <- "sum"
CUpr15$Year <- as.character(CUpr15$Year)
CUpr15_final <- bind_rows(CUpr15, CUpr15_sampled)
CUpr15_final <- CUpr15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBsg15_remove <- CBsg15$Species[(CBsg15$Species %in% CBsg15_10$Species)]
CBsg15_remove <- factor(CBsg15_remove)
# check for actual overlap:
CBsg15_remove
CBsg15_10 <- CBsg15_10 %>%
  filter(!(Species %in% CBsg15_remove))
names(CBsg15_10)[names(CBsg15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBsg15_10_pool <- rep(CBsg15_10$Species, CBsg15_10$sum)
CBsg15_10_pool <- factor(CBsg15_10_pool)
set.seed(933)
CBsg15_sampled <- as.data.frame(table((sample(CBsg15_10_pool, sum(CBsg15_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBsg15_10$sum), length(CBsg15_10_pool))))))
CBsg15_sampled$Site.Name <- "Carrie Bow Seagrass"
CBsg15_sampled$Habitat <- "Seagrass"  
CBsg15_sampled$Year <- "2015"
names(CBsg15_sampled)[names(CBsg15_sampled)=="Var1"] <- "Species"
names(CBsg15_sampled)[names(CBsg15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBsg15)[names(CBsg15)=="sum(Total)"] <- "sum"
CBsg15$Year <- as.character(CBsg15$Year)
CBsg15_final <- bind_rows(CBsg15, CBsg15_sampled)
CBsg15_final <- CBsg15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CUsa15_remove <- CUsa15$Species[(CUsa15$Species %in% CUsa15_10$Species)]
CUsa15_remove <- factor(CUsa15_remove)
# check for actual overlap:
CUsa15_remove
CUsa15_10 <- CUsa15_10 %>%
  filter(!(Species %in% CUsa15_remove))
names(CUsa15_10)[names(CUsa15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUsa15_10_pool <- rep(CUsa15_10$Species, CUsa15_10$sum)
CUsa15_10_pool <- factor(CUsa15_10_pool)
set.seed(572)
CUsa15_sampled <- as.data.frame(table((sample(CUsa15_10_pool, sum(CUsa15_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUsa15_10$sum), length(CUsa15_10_pool))))))
CUsa15_sampled$Site.Name <- "Curlew Sand"
CUsa15_sampled$Habitat <- "Sand"  
CUsa15_sampled$Year <- "2015"
names(CUsa15_sampled)[names(CUsa15_sampled)=="Var1"] <- "Species"
names(CUsa15_sampled)[names(CUsa15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUsa15)[names(CUsa15)=="sum(Total)"] <- "sum"
CUsa15$Year <- as.character(CUsa15$Year)
CUsa15_final <- bind_rows(CUsa15, CUsa15_sampled)
CUsa15_final <- CUsa15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CUsg15_remove <- CUsg15$Species[(CUsg15$Species %in% CUsg15_10$Species)]
CUsg15_remove <- factor(CUsg15_remove)
# check for actual overlap:
CUsg15_remove
CUsg15_10 <- CUsg15_10 %>%
  filter(!(Species %in% CUsg15_remove))
names(CUsg15_10)[names(CUsg15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUsg15_10_pool <- rep(CUsg15_10$Species, CUsg15_10$sum)
CUsg15_10_pool <- factor(CUsg15_10_pool)
set.seed(572)
CUsg15_sampled <- as.data.frame(table((sample(CUsg15_10_pool, sum(CUsg15_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUsg15_10$sum), length(CUsg15_10_pool))))))
CUsg15_sampled$Site.Name <- "Curlew Seagrass"
CUsg15_sampled$Habitat <- "Seagrass"  
CUsg15_sampled$Year <- "2015"
names(CUsg15_sampled)[names(CUsg15_sampled)=="Var1"] <- "Species"
names(CUsg15_sampled)[names(CUsg15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUsg15)[names(CUsg15)=="sum(Total)"] <- "sum"
CUsg15$Year <- as.character(CUsg15$Year)
CUsg15_final <- bind_rows(CUsg15, CUsg15_sampled)
CUsg15_final <- CUsg15_final %>%
  filter(sum != 0)



# remove M1/M2 species from M1
HRpr15_remove <- HRpr15$Species[(HRpr15$Species %in% HRpr15_10$Species)]
HRpr15_remove <- factor(HRpr15_remove)
# check for actual overlap:
HRpr15_remove
HRpr15_10 <- HRpr15_10 %>%
  filter(!(Species %in% HRpr15_remove))
names(HRpr15_10)[names(HRpr15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
HRpr15_10_pool <- rep(HRpr15_10$Species, HRpr15_10$sum)
HRpr15_10_pool <- factor(HRpr15_10_pool)
set.seed(891)
HRpr15_sampled <- as.data.frame(table((sample(HRpr15_10_pool, sum(HRpr15_10$sum)/5, replace = FALSE, prob = rep(1/sum(HRpr15_10$sum), length(HRpr15_10_pool))))))
HRpr15_sampled$Site.Name <- "House Reef"
HRpr15_sampled$Habitat <- "Patch Reef"  
HRpr15_sampled$Year <- "2015"
names(HRpr15_sampled)[names(HRpr15_sampled)=="Var1"] <- "Species"
names(HRpr15_sampled)[names(HRpr15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(HRpr15)[names(HRpr15)=="sum(Total)"] <- "sum"
HRpr15$Year <- as.character(HRpr15$Year)
HRpr15_final <- bind_rows(HRpr15, HRpr15_sampled)
HRpr15_final <- HRpr15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
SRfr15_remove <- SRfr15$Species[(SRfr15$Species %in% SRfr15_10$Species)]
SRfr15_remove <- factor(SRfr15_remove)
# check for actual overlap:
SRfr15_remove
SRfr15_10 <- SRfr15_10 %>%
  filter(!(Species %in% SRfr15_remove))
names(SRfr15_10)[names(SRfr15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
SRfr15_10_pool <- rep(SRfr15_10$Species, SRfr15_10$sum)
SRfr15_10_pool <- factor(SRfr15_10_pool)
set.seed(418)
SRfr15_sampled <- as.data.frame(table((sample(SRfr15_10_pool, sum(SRfr15_10$sum)/5, replace = FALSE, prob = rep(1/sum(SRfr15_10$sum), length(SRfr15_10_pool))))))
SRfr15_sampled$Site.Name <- "South Reef"
SRfr15_sampled$Habitat <- "Fore Reef"  
SRfr15_sampled$Year <- "2015"
names(SRfr15_sampled)[names(SRfr15_sampled)=="Var1"] <- "Species"
names(SRfr15_sampled)[names(SRfr15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(SRfr15)[names(SRfr15)=="sum(Total)"] <- "sum"
SRfr15$Year <- as.character(SRfr15$Year)
SRfr15_final <- bind_rows(SRfr15, SRfr15_sampled)
SRfr15_final <- SRfr15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
TBfr15_remove <- TBfr15$Species[(TBfr15$Species %in% TBfr15_10$Species)]
TBfr15_remove <- factor(TBfr15_remove)
# check for actual overlap:
TBfr15_remove
TBfr15_10 <- TBfr15_10 %>%
  filter(!(Species %in% TBfr15_remove))
names(TBfr15_10)[names(TBfr15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBfr15_10_pool <- rep(TBfr15_10$Species, TBfr15_10$sum)
TBfr15_10_pool <- factor(TBfr15_10_pool)
set.seed(573)
TBfr15_sampled <- as.data.frame(table((sample(TBfr15_10_pool, sum(TBfr15_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBfr15_10$sum), length(TBfr15_10_pool))))))
TBfr15_sampled$Site.Name <- "Tobacco Reef"
TBfr15_sampled$Habitat <- "Fore Reef"  
TBfr15_sampled$Year <- "2015"
names(TBfr15_sampled)[names(TBfr15_sampled)=="Var1"] <- "Species"
names(TBfr15_sampled)[names(TBfr15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBfr15)[names(TBfr15)=="sum(Total)"] <- "sum"
TBfr15$Year <- as.character(TBfr15$Year)
TBfr15_final <- bind_rows(TBfr15, TBfr15_sampled)
TBfr15_final <- TBfr15_final %>%
  filter(sum != 0)



# remove M1/M2 species from M1
TBma15_remove <- TBma15$Species[(TBma15$Species %in% TBma15_10$Species)]
TBma15_remove <- factor(TBma15_remove)
# check for actual overlap:
TBma15_remove
TBma15_10 <- TBma15_10 %>%
  filter(!(Species %in% TBma15_remove))
names(TBma15_10)[names(TBma15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBma15_10_pool <- rep(TBma15_10$Species, TBma15_10$sum)
TBma15_10_pool <- factor(TBma15_10_pool)
set.seed(663)
TBma15_sampled <- as.data.frame(table((sample(TBma15_10_pool, sum(TBma15_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBma15_10$sum), length(TBma15_10_pool))))))
TBma15_sampled$Site.Name <- "Tobacco Mangrove"
TBma15_sampled$Habitat <- "Mangrove"  
TBma15_sampled$Year <- "2015"
names(TBma15_sampled)[names(TBma15_sampled)=="Var1"] <- "Species"
names(TBma15_sampled)[names(TBma15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBma15)[names(TBma15)=="sum(Total)"] <- "sum"
TBma15$Year <- as.character(TBma15$Year)
TBma15_final <- bind_rows(TBma15, TBma15_sampled)
TBma15_final <- TBma15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
TBsa15_remove <- TBsa15$Species[(TBsa15$Species %in% TBsa15_10$Species)]
TBsa15_remove <- factor(TBsa15_remove)
# check for actual overlap:
TBsa15_remove
TBsa15_10 <- TBsa15_10 %>%
  filter(!(Species %in% TBsa15_remove))
names(TBsa15_10)[names(TBsa15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBsa15_10_pool <- rep(TBsa15_10$Species, TBsa15_10$sum)
TBsa15_10_pool <- factor(TBsa15_10_pool)
set.seed(149)
TBsa15_sampled <- as.data.frame(table((sample(TBsa15_10_pool, sum(TBsa15_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBsa15_10$sum), length(TBsa15_10_pool))))))
TBsa15_sampled$Site.Name <- "Tobacco Sand"
TBsa15_sampled$Habitat <- "Sand"  
TBsa15_sampled$Year <- "2015"
names(TBsa15_sampled)[names(TBsa15_sampled)=="Var1"] <- "Species"
names(TBsa15_sampled)[names(TBsa15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBsa15)[names(TBsa15)=="sum(Total)"] <- "sum"
TBsa15$Year <- as.character(TBsa15$Year)
TBsa15_final <- bind_rows(TBsa15, TBsa15_sampled)
TBsa15_final <- TBsa15_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
TBma15_remove <- TBma15$Species[(TBma15$Species %in% TBma15_10$Species)]
TBma15_remove <- factor(TBma15_remove)
# check for actual overlap:
TBma15_remove
TBma15_10 <- TBma15_10 %>%
  filter(!(Species %in% TBma15_remove))
names(TBma15_10)[names(TBma15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBma15_10_pool <- rep(TBma15_10$Species, TBma15_10$sum)
TBma15_10_pool <- factor(TBma15_10_pool)
set.seed(995)
TBma15_sampled <- as.data.frame(table((sample(TBma15_10_pool, sum(TBma15_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBma15_10$sum), length(TBma15_10_pool))))))
TBma15_sampled$Site.Name <- "Tobacco Mangrove"
TBma15_sampled$Habitat <- "Mangrove"  
TBma15_sampled$Year <- "2015"
names(TBma15_sampled)[names(TBma15_sampled)=="Var1"] <- "Species"
names(TBma15_sampled)[names(TBma15_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBma15)[names(TBma15)=="sum(Total)"] <- "sum"
TBma15$Year <- as.character(TBma15$Year)
TBma15_final <- bind_rows(TBma15, TBma15_sampled)
TBma15_final <- TBma15_final %>%
  filter(sum != 0)



# No M2 species
names(TCma15_10)[names(TCma15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TCma15_10_pool <- rep(TCma15_10$Species, TCma15_10$sum)
TCma15_10_pool <- factor(TCma15_10_pool)
set.seed(968)
TCma15_sampled <- as.data.frame(table((sample(TCma15_10_pool, sum(TCma15_10$sum)/5, replace = FALSE, prob = rep(1/sum(TCma15_10$sum), length(TCma15_10_pool))))))
TCma15_sampled$Site.Name <- "Twin Cays Mangrove"
TCma15_sampled$Habitat <- "Mangrove"  
TCma15_sampled$Year <- "2015"
names(TCma15_sampled)[names(TCma15_sampled)=="Var1"] <- "Species"
names(TCma15_sampled)[names(TCma15_sampled)=="Freq"] <- "sum"
# filter zero draws
TCma15_final <- TCma15_sampled %>%
  filter(sum != 0)



# No M2 species
names(TCsg15_10)[names(TCsg15_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TCsg15_10_pool <- rep(TCsg15_10$Species, TCsg15_10$sum)
TCsg15_10_pool <- factor(TCsg15_10_pool)
set.seed(559)
TCsg15_sampled <- as.data.frame(table((sample(TCsg15_10_pool, sum(TCsg15_10$sum)/5, replace = FALSE, prob = rep(1/sum(TCsg15_10$sum), length(TCsg15_10_pool))))))
TCsg15_sampled$Site.Name <- "Twin Cays Seagrass"
TCsg15_sampled$Habitat <- "Seagrass"  
TCsg15_sampled$Year <- "2015"
names(TCsg15_sampled)[names(TCsg15_sampled)=="Var1"] <- "Species"
names(TCsg15_sampled)[names(TCsg15_sampled)=="Freq"] <- "sum"
# filter zero draws
TCsg15_final <- TCsg15_sampled %>%
  filter(sum != 0)

########
# 2016 #
########

# No M2 species
names(BGma16_10)[names(BGma16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
BGma16_10_pool <- rep(BGma16_10$Species, BGma16_10$sum)
BGma16_10_pool <- factor(BGma16_10_pool)
set.seed(396)
BGma16_sampled <- as.data.frame(table((sample(BGma16_10_pool, sum(BGma16_10$sum)/5, replace = FALSE, prob = rep(1/sum(BGma16_10$sum), length(BGma16_10_pool))))))
BGma16_sampled$Site.Name <- "Blueground Mangrove"
BGma16_sampled$Habitat <- "Mangrove"  
BGma16_sampled$Year <- "2016"
names(BGma16_sampled)[names(BGma16_sampled)=="Var1"] <- "Species"
names(BGma16_sampled)[names(BGma16_sampled)=="Freq"] <- "sum"
# filter zero draws
BGma16_final <- BGma16_sampled %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBfr16_remove <- CBfr16$Species[(CBfr16$Species %in% CBfr16_10$Species)]
CBfr16_remove <- factor(CBfr16_remove)
# check for actual overlap:
CBfr16_remove
CBfr16_10 <- CBfr16_10 %>%
  filter(!(Species %in% CBfr16_remove))
names(CBfr16_10)[names(CBfr16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBfr16_10_pool <- rep(CBfr16_10$Species, CBfr16_10$sum)
CBfr16_10_pool <- factor(CBfr16_10_pool)
set.seed(967)
CBfr16_sampled <- as.data.frame(table((sample(CBfr16_10_pool, sum(CBfr16_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBfr16_10$sum), length(CBfr16_10_pool))))))
CBfr16_sampled$Site.Name <- "Carrie Bow Reef"
CBfr16_sampled$Habitat <- "Fore Reef"  
CBfr16_sampled$Year <- "2016"
names(CBfr16_sampled)[names(CBfr16_sampled)=="Var1"] <- "Species"
names(CBfr16_sampled)[names(CBfr16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBfr16)[names(CBfr16)=="sum(Total)"] <- "sum"
CBfr16$Year <- as.character(CBfr16$Year)
CBfr16_final <- bind_rows(CBfr16, CBfr16_sampled)
CBfr16_final <- CBfr16_final %>%
  filter(sum != 0)
# harmonize name
library(plyr)
CBfr16_final$Site.Name <- CBfr16_final$Site.Name %>%
  revalue(c("Carrie Bow Reef" = "Carrie Bow Reef"))
detach(package:plyr)



# remove M1/M2 species from M1
CBpa16_remove <- CBpa16$Species[(CBpa16$Species %in% CBpa16_10$Species)]
CBpa16_remove <- factor(CBpa16_remove)
# check for actual overlap:
CBpa16_remove
CBpa16_10 <- CBpa16_10 %>%
  filter(!(Species %in% CBpa16_remove))
names(CBpa16_10)[names(CBpa16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBpa16_10_pool <- rep(CBpa16_10$Species, CBpa16_10$sum)
CBpa16_10_pool <- factor(CBpa16_10_pool)
set.seed(692)
CBpa16_sampled <- as.data.frame(table((sample(CBpa16_10_pool, sum(CBpa16_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBpa16_10$sum), length(CBpa16_10_pool))))))
CBpa16_sampled$Site.Name <- "Carrie Bow Patch Reef"
CBpa16_sampled$Habitat <- "Patch Reef"  
CBpa16_sampled$Year <- "2016"
names(CBpa16_sampled)[names(CBpa16_sampled)=="Var1"] <- "Species"
names(CBpa16_sampled)[names(CBpa16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBpa16)[names(CBpa16)=="sum(Total)"] <- "sum"
CBpa16$Year <- as.character(CBpa16$Year)
CBpa16_final <- bind_rows(CBpa16, CBpa16_sampled)
CBpa16_final <- CBpa16_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBsa16_remove <- CBsa16$Species[(CBsa16$Species %in% CBsa16_10$Species)]
CBsa16_remove <- factor(CBsa16_remove)
# check for actual overlap:
CBsa16_remove
CBsa16_10 <- CBsa16_10 %>%
  filter(!(Species %in% CBsa16_remove))
names(CBsa16_10)[names(CBsa16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBsa16_10_pool <- rep(CBsa16_10$Species, CBsa16_10$sum)
CBsa16_10_pool <- factor(CBsa16_10_pool)
set.seed(578)
CBsa16_sampled <- as.data.frame(table((sample(CBsa16_10_pool, sum(CBsa16_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBsa16_10$sum), length(CBsa16_10_pool))))))
CBsa16_sampled$Site.Name <- "Carrie Bow Sand"
CBsa16_sampled$Habitat <- "Sand"  
CBsa16_sampled$Year <- "2016"
names(CBsa16_sampled)[names(CBsa16_sampled)=="Var1"] <- "Species"
names(CBsa16_sampled)[names(CBsa16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBsa16)[names(CBsa16)=="sum(Total)"] <- "sum"
CBsa16$Year <- as.character(CBsa16$Year)
CBsa16_final <- bind_rows(CBsa16, CBsa16_sampled)
CBsa16_final <- CBsa16_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBsg16_remove <- CBsg16$Species[(CBsg16$Species %in% CBsg16_10$Species)]
CBsg16_remove <- factor(CBsg16_remove)
# check for actual overlap:
CBsg16_remove
CBsg16_10 <- CBsg16_10 %>%
  filter(!(Species %in% CBsg16_remove))
names(CBsg16_10)[names(CBsg16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBsg16_10_pool <- rep(CBsg16_10$Species, CBsg16_10$sum)
CBsg16_10_pool <- factor(CBsg16_10_pool)
set.seed(262)
CBsg16_sampled <- as.data.frame(table((sample(CBsg16_10_pool, sum(CBsg16_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBsg16_10$sum), length(CBsg16_10_pool))))))
CBsg16_sampled$Site.Name <- "Carrie Bow Seagrass"
CBsg16_sampled$Habitat <- "Seagrass"  
CBsg16_sampled$Year <- "2016"
names(CBsg16_sampled)[names(CBsg16_sampled)=="Var1"] <- "Species"
names(CBsg16_sampled)[names(CBsg16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBsg16)[names(CBsg16)=="sum(Total)"] <- "sum"
CBsg16$Year <- as.character(CBsg16$Year)
CBsg16_final <- bind_rows(CBsg16, CBsg16_sampled)
CBsg16_final <- CBsg16_final %>%
  filter(sum != 0)



# remove M1/M2 species from M1
CUpr16_remove <- CUpr16$Species[(CUpr16$Species %in% CUpr16_10$Species)]
CUpr16_remove <- factor(CUpr16_remove)
# check for actual overlap:
CUpr16_remove
CUpr16_10 <- CUpr16_10 %>%
  filter(!(Species %in% CUpr16_remove))
names(CUpr16_10)[names(CUpr16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUpr16_10_pool <- rep(CUpr16_10$Species, CUpr16_10$sum)
CUpr16_10_pool <- factor(CUpr16_10_pool)
set.seed(983)
CUpr16_sampled <- as.data.frame(table((sample(CUpr16_10_pool, sum(CUpr16_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUpr16_10$sum), length(CUpr16_10_pool))))))
CUpr16_sampled$Site.Name <- "Curlew Patch Reef"
CUpr16_sampled$Habitat <- "Patch Reef"  
CUpr16_sampled$Year <- "2016"
names(CUpr16_sampled)[names(CUpr16_sampled)=="Var1"] <- "Species"
names(CUpr16_sampled)[names(CUpr16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUpr16)[names(CUpr16)=="sum(Total)"] <- "sum"
CUpr16$Year <- as.character(CUpr16$Year)
CUpr16_final <- bind_rows(CUpr16, CUpr16_sampled)
CUpr16_final <- CUpr16_final %>%
  filter(sum != 0)



# remove M1/M2 species from M1
CUsa16_remove <- CUsa16$Species[(CUsa16$Species %in% CUsa16_10$Species)]
CUsa16_remove <- factor(CUsa16_remove)
# check for actual overlap:
CUsa16_remove
CUsa16_10 <- CUsa16_10 %>%
  filter(!(Species %in% CUsa16_remove))
names(CUsa16_10)[names(CUsa16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUsa16_10_pool <- rep(CUsa16_10$Species, CUsa16_10$sum)
CUsa16_10_pool <- factor(CUsa16_10_pool)
set.seed(652)
CUsa16_sampled <- as.data.frame(table((sample(CUsa16_10_pool, sum(CUsa16_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUsa16_10$sum), length(CUsa16_10_pool))))))
CUsa16_sampled$Site.Name <- "Curlew Sand"
CUsa16_sampled$Habitat <- "Sand"  
CUsa16_sampled$Year <- "2016"
names(CUsa16_sampled)[names(CUsa16_sampled)=="Var1"] <- "Species"
names(CUsa16_sampled)[names(CUsa16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUsa16)[names(CUsa16)=="sum(Total)"] <- "sum"
CUsa16$Year <- as.character(CUsa16$Year)
CUsa16_final <- bind_rows(CUsa16, CUsa16_sampled)
CUsa16_final <- CUsa16_final %>%
  filter(sum != 0)




# remove M1/M2 species from M1
CUsg16_remove <- CUsg16$Species[(CUsg16$Species %in% CUsg16_10$Species)]
CUsg16_remove <- factor(CUsg16_remove)
# check for actual overlap:
CUsg16_remove
CUsg16_10 <- CUsg16_10 %>%
  filter(!(Species %in% CUsg16_remove))
names(CUsg16_10)[names(CUsg16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUsg16_10_pool <- rep(CUsg16_10$Species, CUsg16_10$sum)
CUsg16_10_pool <- factor(CUsg16_10_pool)
set.seed(369)
CUsg16_sampled <- as.data.frame(table((sample(CUsg16_10_pool, sum(CUsg16_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUsg16_10$sum), length(CUsg16_10_pool))))))
CUsg16_sampled$Site.Name <- "Curlew Seagrass"
CUsg16_sampled$Habitat <- "Seagrass"  
CUsg16_sampled$Year <- "2016"
names(CUsg16_sampled)[names(CUsg16_sampled)=="Var1"] <- "Species"
names(CUsg16_sampled)[names(CUsg16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUsg16)[names(CUsg16)=="sum(Total)"] <- "sum"
CUsg16$Year <- as.character(CUsg16$Year)
CUsg16_final <- bind_rows(CUsg16, CUsg16_sampled)
CUsg16_final <- CUsg16_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
HRpr16_remove <- HRpr16$Species[(HRpr16$Species %in% HRpr16_10$Species)]
HRpr16_remove <- factor(HRpr16_remove)
# check for actual overlap:
HRpr16_remove
HRpr16_10 <- HRpr16_10 %>%
  filter(!(Species %in% HRpr16_remove))
names(HRpr16_10)[names(HRpr16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
HRpr16_10_pool <- rep(HRpr16_10$Species, HRpr16_10$sum)
HRpr16_10_pool <- factor(HRpr16_10_pool)
set.seed(611)
HRpr16_sampled <- as.data.frame(table((sample(HRpr16_10_pool, sum(HRpr16_10$sum)/5, replace = FALSE, prob = rep(1/sum(HRpr16_10$sum), length(HRpr16_10_pool))))))
HRpr16_sampled$Site.Name <- "House Reef"
HRpr16_sampled$Habitat <- "Patch Reef"  
HRpr16_sampled$Year <- "2016"
names(HRpr16_sampled)[names(HRpr16_sampled)=="Var1"] <- "Species"
names(HRpr16_sampled)[names(HRpr16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(HRpr16)[names(HRpr16)=="sum(Total)"] <- "sum"
HRpr16$Year <- as.character(HRpr16$Year)
HRpr16_final <- bind_rows(HRpr16, HRpr16_sampled)
HRpr16_final <- HRpr16_final %>%
  filter(sum != 0)



# remove M1/M2 species from M1
SRfr16_remove <- SRfr16$Species[(SRfr16$Species %in% SRfr16_10$Species)]
SRfr16_remove <- factor(SRfr16_remove)
# check for actual overlap:
SRfr16_remove
SRfr16_10 <- SRfr16_10 %>%
  filter(!(Species %in% SRfr16_remove))
names(SRfr16_10)[names(SRfr16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
SRfr16_10_pool <- rep(SRfr16_10$Species, SRfr16_10$sum)
SRfr16_10_pool <- factor(SRfr16_10_pool)
set.seed(184)
SRfr16_sampled <- as.data.frame(table((sample(SRfr16_10_pool, sum(SRfr16_10$sum)/5, replace = FALSE, prob = rep(1/sum(SRfr16_10$sum), length(SRfr16_10_pool))))))
SRfr16_sampled$Site.Name <- "South Reef"
SRfr16_sampled$Habitat <- "Fore Reef"  
SRfr16_sampled$Year <- "2016"
names(SRfr16_sampled)[names(SRfr16_sampled)=="Var1"] <- "Species"
names(SRfr16_sampled)[names(SRfr16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(SRfr16)[names(SRfr16)=="sum(Total)"] <- "sum"
SRfr16$Year <- as.character(SRfr16$Year)
SRfr16_final <- bind_rows(SRfr16, SRfr16_sampled)
SRfr16_final <- SRfr16_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
TBfr16_remove <- TBfr16$Species[(TBfr16$Species %in% TBfr16_10$Species)]
TBfr16_remove <- factor(TBfr16_remove)
# check for actual overlap:
TBfr16_remove
TBfr16_10 <- TBfr16_10 %>%
  filter(!(Species %in% TBfr16_remove))
names(TBfr16_10)[names(TBfr16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBfr16_10_pool <- rep(TBfr16_10$Species, TBfr16_10$sum)
TBfr16_10_pool <- factor(TBfr16_10_pool)
set.seed(290)
TBfr16_sampled <- as.data.frame(table((sample(TBfr16_10_pool, sum(TBfr16_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBfr16_10$sum), length(TBfr16_10_pool))))))
TBfr16_sampled$Site.Name <- "Tobacco Reef"
TBfr16_sampled$Habitat <- "Fore Reef"  
TBfr16_sampled$Year <- "2016"
names(TBfr16_sampled)[names(TBfr16_sampled)=="Var1"] <- "Species"
names(TBfr16_sampled)[names(TBfr16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBfr16)[names(TBfr16)=="sum(Total)"] <- "sum"
TBfr16$Year <- as.character(TBfr16$Year)
TBfr16_final <- bind_rows(TBfr16, TBfr16_sampled)
TBfr16_final <- TBfr16_final %>%
  filter(sum != 0)



# remove M1/M2 species from M1
TBma16_remove <- TBma16$Species[(TBma16$Species %in% TBma16_10$Species)]
TBma16_remove <- factor(TBma16_remove)
# check for actual overlap:
TBma16_remove
TBma16_10 <- TBma16_10 %>%
  filter(!(Species %in% TBma16_remove))
names(TBma16_10)[names(TBma16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBma16_10_pool <- rep(TBma16_10$Species, TBma16_10$sum)
TBma16_10_pool <- factor(TBma16_10_pool)
set.seed(965)
TBma16_sampled <- as.data.frame(table((sample(TBma16_10_pool, sum(TBma16_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBma16_10$sum), length(TBma16_10_pool))))))
TBma16_sampled$Site.Name <- "Tobacco Mangrove"
TBma16_sampled$Habitat <- "Mangrove"  
TBma16_sampled$Year <- "2016"
names(TBma16_sampled)[names(TBma16_sampled)=="Var1"] <- "Species"
names(TBma16_sampled)[names(TBma16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBma16)[names(TBma16)=="sum(Total)"] <- "sum"
TBma16$Year <- as.character(TBma16$Year)
TBma16_final <- bind_rows(TBma16, TBma16_sampled)
TBma16_final <- TBma16_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
TBsa16_remove <- TBsa16$Species[(TBsa16$Species %in% TBsa16_10$Species)]
TBsa16_remove <- factor(TBsa16_remove)
# check for actual overlap:
TBsa16_remove
TBsa16_10 <- TBsa16_10 %>%
  filter(!(Species %in% TBsa16_remove))
names(TBsa16_10)[names(TBsa16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBsa16_10_pool <- rep(TBsa16_10$Species, TBsa16_10$sum)
TBsa16_10_pool <- factor(TBsa16_10_pool)
set.seed(871)
TBsa16_sampled <- as.data.frame(table((sample(TBsa16_10_pool, sum(TBsa16_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBsa16_10$sum), length(TBsa16_10_pool))))))
TBsa16_sampled$Site.Name <- "Tobacco Sand"
TBsa16_sampled$Habitat <- "Sand"  
TBsa16_sampled$Year <- "2016"
names(TBsa16_sampled)[names(TBsa16_sampled)=="Var1"] <- "Species"
names(TBsa16_sampled)[names(TBsa16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBsa16)[names(TBsa16)=="sum(Total)"] <- "sum"
TBsa16$Year <- as.character(TBsa16$Year)
TBsa16_final <- bind_rows(TBsa16, TBsa16_sampled)
TBsa16_final <- TBsa16_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBpa16_remove <- CBpa16$Species[(CBpa16$Species %in% CBpa16_10$Species)]
CBpa16_remove <- factor(CBpa16_remove)
# check for actual overlap:
CBpa16_remove
CBpa16_10 <- CBpa16_10 %>%
  filter(!(Species %in% CBpa16_remove))
names(CBpa16_10)[names(CBpa16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBpa16_10_pool <- rep(CBpa16_10$Species, CBpa16_10$sum)
CBpa16_10_pool <- factor(CBpa16_10_pool)
set.seed(451)
CBpa16_sampled <- as.data.frame(table((sample(CBpa16_10_pool, sum(CBpa16_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBpa16_10$sum), length(CBpa16_10_pool))))))
CBpa16_sampled$Site.Name <- "Carrie Bow Patch Reef"
CBpa16_sampled$Habitat <- "Patch Reef"  
CBpa16_sampled$Year <- "2016"
names(CBpa16_sampled)[names(CBpa16_sampled)=="Var1"] <- "Species"
names(CBpa16_sampled)[names(CBpa16_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBpa16)[names(CBpa16)=="sum(Total)"] <- "sum"
CBpa16$Year <- as.character(CBpa16$Year)
CBpa16_final <- bind_rows(CBpa16, CBpa16_sampled)
CBpa16_final <- CBpa16_final %>%
  filter(sum != 0)


# No M2 species
names(TCma16_10)[names(TCma16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TCma16_10_pool <- rep(TCma16_10$Species, TCma16_10$sum)
TCma16_10_pool <- factor(TCma16_10_pool)
set.seed(521)
TCma16_sampled <- as.data.frame(table((sample(TCma16_10_pool, sum(TCma16_10$sum)/5, replace = FALSE, prob = rep(1/sum(TCma16_10$sum), length(TCma16_10_pool))))))
TCma16_sampled$Site.Name <- "Twin Cays Mangrove"
TCma16_sampled$Habitat <- "Mangrove"  
TCma16_sampled$Year <- "2016"
names(TCma16_sampled)[names(TCma16_sampled)=="Var1"] <- "Species"
names(TCma16_sampled)[names(TCma16_sampled)=="Freq"] <- "sum"
# filter zero draws
TCma16_final <- TCma16_sampled %>%
  filter(sum != 0)


# No M2 species
names(TCsg16_10)[names(TCsg16_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TCsg16_10_pool <- rep(TCsg16_10$Species, TCsg16_10$sum)
TCsg16_10_pool <- factor(TCsg16_10_pool)
set.seed(682)
TCsg16_sampled <- as.data.frame(table((sample(TCsg16_10_pool, sum(TCsg16_10$sum)/5, replace = FALSE, prob = rep(1/sum(TCsg16_10$sum), length(TCsg16_10_pool))))))
TCsg16_sampled$Site.Name <- "Twin Cays Seagrass"
TCsg16_sampled$Habitat <- "Seagrass"  
TCsg16_sampled$Year <- "2016"
names(TCsg16_sampled)[names(TCsg16_sampled)=="Var1"] <- "Species"
names(TCsg16_sampled)[names(TCsg16_sampled)=="Freq"] <- "sum"
# filter zero draws
TCsg16_final <- TCsg16_sampled %>%
  filter(sum != 0)

########
# 2017 #
########

# no M2 species
names(BGma17_10)[names(BGma17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
BGma17_10_pool <- rep(BGma17_10$Species, BGma17_10$sum)
BGma17_10_pool <- factor(BGma17_10_pool)
set.seed(396)
BGma17_sampled <- as.data.frame(table((sample(BGma17_10_pool, sum(BGma17_10$sum)/5, replace = FALSE, prob = rep(1/sum(BGma17_10$sum), length(BGma17_10_pool))))))
BGma17_sampled$Site.Name <- "Blueground Mangrove"
BGma17_sampled$Habitat <- "Mangrove"  
BGma17_sampled$Year <- "2017"
names(BGma17_sampled)[names(BGma17_sampled)=="Var1"] <- "Species"
names(BGma17_sampled)[names(BGma17_sampled)=="Freq"] <- "sum"
# filter zero draws
BGma17_final <- BGma17_sampled %>%
  filter(sum != 0)


# no M2 species
names(CBfr17_10)[names(CBfr17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBfr17_10_pool <- rep(CBfr17_10$Species, CBfr17_10$sum)
CBfr17_10_pool <- factor(CBfr17_10_pool)
set.seed(967)
CBfr17_sampled <- as.data.frame(table((sample(CBfr17_10_pool, sum(CBfr17_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBfr17_10$sum), length(CBfr17_10_pool))))))
CBfr17_sampled$Site.Name <- "Carrie Bow Reef"
CBfr17_sampled$Habitat <- "Fore Reef"  
CBfr17_sampled$Year <- "2017"
names(CBfr17_sampled)[names(CBfr17_sampled)=="Var1"] <- "Species"
names(CBfr17_sampled)[names(CBfr17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBfr17)[names(CBfr17)=="sum(Total)"] <- "sum"
CBfr17$Year <- as.character(CBfr17$Year)
CBfr17_final <- bind_rows(CBfr17, CBfr17_sampled)
CBfr17_final <- CBfr17_final %>%
  filter(sum != 0)


# no M2 species
names(CBpa17_10)[names(CBpa17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBpa17_10_pool <- rep(CBpa17_10$Species, CBpa17_10$sum)
CBpa17_10_pool <- factor(CBpa17_10_pool)
set.seed(692)
CBpa17_sampled <- as.data.frame(table((sample(CBpa17_10_pool, sum(CBpa17_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBpa17_10$sum), length(CBpa17_10_pool))))))
CBpa17_sampled$Site.Name <- "Carrie Bow Patch Reef"
CBpa17_sampled$Habitat <- "Patch Reef"  
CBpa17_sampled$Year <- "2017"
names(CBpa17_sampled)[names(CBpa17_sampled)=="Var1"] <- "Species"
names(CBpa17_sampled)[names(CBpa17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBpa17)[names(CBpa17)=="sum(Total)"] <- "sum"
CBpa17$Year <- as.character(CBpa17$Year)
CBpa17_final <- bind_rows(CBpa17, CBpa17_sampled)
CBpa17_final <- CBpa17_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBsa17_remove <- CBsa17$Species[(CBsa17$Species %in% CBsa17_10$Species)]
CBsa17_remove <- factor(CBsa17_remove)
# check for actual overlap:
CBsa17_remove
CBsa17_10 <- CBsa17_10 %>%
  filter(!(Species %in% CBsa17_remove))
names(CBsa17_10)[names(CBsa17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBsa17_10_pool <- rep(CBsa17_10$Species, CBsa17_10$sum)
CBsa17_10_pool <- factor(CBsa17_10_pool)
set.seed(578)
CBsa17_sampled <- as.data.frame(table((sample(CBsa17_10_pool, sum(CBsa17_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBsa17_10$sum), length(CBsa17_10_pool))))))
CBsa17_sampled$Site.Name <- "Carrie Bow Sand"
CBsa17_sampled$Habitat <- "Sand"  
CBsa17_sampled$Year <- "2017"
names(CBsa17_sampled)[names(CBsa17_sampled)=="Var1"] <- "Species"
names(CBsa17_sampled)[names(CBsa17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBsa17)[names(CBsa17)=="sum(Total)"] <- "sum"
CBsa17$Year <- as.character(CBsa17$Year)
CBsa17_final <- bind_rows(CBsa17, CBsa17_sampled)
CBsa17_final <- CBsa17_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
CBsg17_remove <- CBsg17$Species[(CBsg17$Species %in% CBsg17_10$Species)]
CBsg17_remove <- factor(CBsg17_remove)
# check for actual overlap:
CBsg17_remove
CBsg17_10 <- CBsg17_10 %>%
  filter(!(Species %in% CBsg17_remove))
names(CBsg17_10)[names(CBsg17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CBsg17_10_pool <- rep(CBsg17_10$Species, CBsg17_10$sum)
CBsg17_10_pool <- factor(CBsg17_10_pool)
set.seed(262)
CBsg17_sampled <- as.data.frame(table((sample(CBsg17_10_pool, sum(CBsg17_10$sum)/5, replace = FALSE, prob = rep(1/sum(CBsg17_10$sum), length(CBsg17_10_pool))))))
CBsg17_sampled$Site.Name <- "Carrie Bow Seagrass"
CBsg17_sampled$Habitat <- "Seagrass"  
CBsg17_sampled$Year <- "2017"
names(CBsg17_sampled)[names(CBsg17_sampled)=="Var1"] <- "Species"
names(CBsg17_sampled)[names(CBsg17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CBsg17)[names(CBsg17)=="sum(Total)"] <- "sum"
CBsg17$Year <- as.character(CBsg17$Year)
CBsg17_final <- bind_rows(CBsg17, CBsg17_sampled)
CBsg17_final <- CBsg17_final %>%
  filter(sum != 0)



# no M2 species
names(CUpr17_10)[names(CUpr17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUpr17_10_pool <- rep(CUpr17_10$Species, CUpr17_10$sum)
CUpr17_10_pool <- factor(CUpr17_10_pool)
set.seed(983)
CUpr17_sampled <- as.data.frame(table((sample(CUpr17_10_pool, sum(CUpr17_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUpr17_10$sum), length(CUpr17_10_pool))))))
CUpr17_sampled$Site.Name <- "Curlew Patch Reef"
CUpr17_sampled$Habitat <- "Patch Reef"  
CUpr17_sampled$Year <- "2017"
names(CUpr17_sampled)[names(CUpr17_sampled)=="Var1"] <- "Species"
names(CUpr17_sampled)[names(CUpr17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUpr17)[names(CUpr17)=="sum(Total)"] <- "sum"
CUpr17$Year <- as.character(CUpr17$Year)
CUpr17_final <- bind_rows(CUpr17, CUpr17_sampled)
CUpr17_final <- CUpr17_final %>%
  filter(sum != 0)



# no M2 species
names(CUsa17_10)[names(CUsa17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUsa17_10_pool <- rep(CUsa17_10$Species, CUsa17_10$sum)
CUsa17_10_pool <- factor(CUsa17_10_pool)
set.seed(652)
CUsa17_sampled <- as.data.frame(table((sample(CUsa17_10_pool, sum(CUsa17_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUsa17_10$sum), length(CUsa17_10_pool))))))
CUsa17_sampled$Site.Name <- "Curlew Sand"
CUsa17_sampled$Habitat <- "Sand"  
CUsa17_sampled$Year <- "2017"
names(CUsa17_sampled)[names(CUsa17_sampled)=="Var1"] <- "Species"
names(CUsa17_sampled)[names(CUsa17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUsa17)[names(CUsa17)=="sum(Total)"] <- "sum"
CUsa17$Year <- as.character(CUsa17$Year)
CUsa17_final <- bind_rows(CUsa17, CUsa17_sampled)
CUsa17_final <- CUsa17_final %>%
  filter(sum != 0)




# no M2 species
names(CUsg17_10)[names(CUsg17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
CUsg17_10_pool <- rep(CUsg17_10$Species, CUsg17_10$sum)
CUsg17_10_pool <- factor(CUsg17_10_pool)
set.seed(369)
CUsg17_sampled <- as.data.frame(table((sample(CUsg17_10_pool, sum(CUsg17_10$sum)/5, replace = FALSE, prob = rep(1/sum(CUsg17_10$sum), length(CUsg17_10_pool))))))
CUsg17_sampled$Site.Name <- "Curlew Seagrass"
CUsg17_sampled$Habitat <- "Seagrass"  
CUsg17_sampled$Year <- "2017"
names(CUsg17_sampled)[names(CUsg17_sampled)=="Var1"] <- "Species"
names(CUsg17_sampled)[names(CUsg17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(CUsg17)[names(CUsg17)=="sum(Total)"] <- "sum"
CUsg17$Year <- as.character(CUsg17$Year)
CUsg17_final <- bind_rows(CUsg17, CUsg17_sampled)
CUsg17_final <- CUsg17_final %>%
  filter(sum != 0)


# no M2 species
names(HRpr17_10)[names(HRpr17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
HRpr17_10_pool <- rep(HRpr17_10$Species, HRpr17_10$sum)
HRpr17_10_pool <- factor(HRpr17_10_pool)
set.seed(611)
HRpr17_sampled <- as.data.frame(table((sample(HRpr17_10_pool, sum(HRpr17_10$sum)/5, replace = FALSE, prob = rep(1/sum(HRpr17_10$sum), length(HRpr17_10_pool))))))
HRpr17_sampled$Site.Name <- "House Reef"
HRpr17_sampled$Habitat <- "Patch Reef"  
HRpr17_sampled$Year <- "2017"
names(HRpr17_sampled)[names(HRpr17_sampled)=="Var1"] <- "Species"
names(HRpr17_sampled)[names(HRpr17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(HRpr17)[names(HRpr17)=="sum(Total)"] <- "sum"
HRpr17$Year <- as.character(HRpr17$Year)
HRpr17_final <- bind_rows(HRpr17, HRpr17_sampled)
HRpr17_final <- HRpr17_final %>%
  filter(sum != 0)



# remove M1/M2 species from M1
SRfr17_remove <- SRfr17$Species[(SRfr17$Species %in% SRfr17_10$Species)]
SRfr17_remove <- factor(SRfr17_remove)
# check for actual overlap:
SRfr17_remove
SRfr17_10 <- SRfr17_10 %>%
  filter(!(Species %in% SRfr17_remove))
names(SRfr17_10)[names(SRfr17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
SRfr17_10_pool <- rep(SRfr17_10$Species, SRfr17_10$sum)
SRfr17_10_pool <- factor(SRfr17_10_pool)
set.seed(184)
SRfr17_sampled <- as.data.frame(table((sample(SRfr17_10_pool, sum(SRfr17_10$sum)/5, replace = FALSE, prob = rep(1/sum(SRfr17_10$sum), length(SRfr17_10_pool))))))
SRfr17_sampled$Site.Name <- "South Reef"
SRfr17_sampled$Habitat <- "Fore Reef"  
SRfr17_sampled$Year <- "2017"
names(SRfr17_sampled)[names(SRfr17_sampled)=="Var1"] <- "Species"
names(SRfr17_sampled)[names(SRfr17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(SRfr17)[names(SRfr17)=="sum(Total)"] <- "sum"
SRfr17$Year <- as.character(SRfr17$Year)
SRfr17_final <- bind_rows(SRfr17, SRfr17_sampled)
SRfr17_final <- SRfr17_final %>%
  filter(sum != 0)


# remove M1/M2 species from M1
TBfr17_remove <- TBfr17$Species[(TBfr17$Species %in% TBfr17_10$Species)]
TBfr17_remove <- factor(TBfr17_remove)
# check for actual overlap:
TBfr17_remove
TBfr17_10 <- TBfr17_10 %>%
  filter(!(Species %in% TBfr17_remove))
names(TBfr17_10)[names(TBfr17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBfr17_10_pool <- rep(TBfr17_10$Species, TBfr17_10$sum)
TBfr17_10_pool <- factor(TBfr17_10_pool)
set.seed(290)
TBfr17_sampled <- as.data.frame(table((sample(TBfr17_10_pool, sum(TBfr17_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBfr17_10$sum), length(TBfr17_10_pool))))))
TBfr17_sampled$Site.Name <- "Tobacco Reef"
TBfr17_sampled$Habitat <- "Fore Reef"  
TBfr17_sampled$Year <- "2017"
names(TBfr17_sampled)[names(TBfr17_sampled)=="Var1"] <- "Species"
names(TBfr17_sampled)[names(TBfr17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBfr17)[names(TBfr17)=="sum(Total)"] <- "sum"
TBfr17$Year <- as.character(TBfr17$Year)
TBfr17_final <- bind_rows(TBfr17, TBfr17_sampled)
TBfr17_final <- TBfr17_final %>%
  filter(sum != 0)



# remove M1/M2 species from M1
TBma17_remove <- TBma17$Species[(TBma17$Species %in% TBma17_10$Species)]
TBma17_remove <- factor(TBma17_remove)
# check for actual overlap:
TBma17_remove
TBma17_10 <- TBma17_10 %>%
  filter(!(Species %in% TBma17_remove))
names(TBma17_10)[names(TBma17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBma17_10_pool <- rep(TBma17_10$Species, TBma17_10$sum)
TBma17_10_pool <- factor(TBma17_10_pool)
set.seed(965)
TBma17_sampled <- as.data.frame(table((sample(TBma17_10_pool, sum(TBma17_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBma17_10$sum), length(TBma17_10_pool))))))
TBma17_sampled$Site.Name <- "Tobacco Mangrove"
TBma17_sampled$Habitat <- "Mangrove"  
TBma17_sampled$Year <- "2017"
names(TBma17_sampled)[names(TBma17_sampled)=="Var1"] <- "Species"
names(TBma17_sampled)[names(TBma17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBma17)[names(TBma17)=="sum(Total)"] <- "sum"
TBma17$Year <- as.character(TBma17$Year)
TBma17_final <- bind_rows(TBma17, TBma17_sampled)
TBma17_final <- TBma17_final %>%
  filter(sum != 0)


# no M2 species
names(TBsa17_10)[names(TBsa17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBsa17_10_pool <- rep(TBsa17_10$Species, TBsa17_10$sum)
TBsa17_10_pool <- factor(TBsa17_10_pool)
set.seed(871)
TBsa17_sampled <- as.data.frame(table((sample(TBsa17_10_pool, sum(TBsa17_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBsa17_10$sum), length(TBsa17_10_pool))))))
TBsa17_sampled$Site.Name <- "Tobacco Sand"
TBsa17_sampled$Habitat <- "Sand"  
TBsa17_sampled$Year <- "2017"
names(TBsa17_sampled)[names(TBsa17_sampled)=="Var1"] <- "Species"
names(TBsa17_sampled)[names(TBsa17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBsa17)[names(TBsa17)=="sum(Total)"] <- "sum"
TBsa17$Year <- as.character(TBsa17$Year)
TBsa17_final <- bind_rows(TBsa17, TBsa17_sampled)
TBsa17_final <- TBsa17_final %>%
  filter(sum != 0)


# no M2 species
names(TBsg17_10)[names(TBsg17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TBsg17_10_pool <- rep(TBsg17_10$Species, TBsg17_10$sum)
TBsg17_10_pool <- factor(TBsg17_10_pool)
set.seed(451)
TBsg17_sampled <- as.data.frame(table((sample(TBsg17_10_pool, sum(TBsg17_10$sum)/5, replace = FALSE, prob = rep(1/sum(TBsg17_10$sum), length(TBsg17_10_pool))))))
TBsg17_sampled$Site.Name <- "Tobacco Seagrass"
TBsg17_sampled$Habitat <- "Seagrass"  
TBsg17_sampled$Year <- "2017"
names(TBsg17_sampled)[names(TBsg17_sampled)=="Var1"] <- "Species"
names(TBsg17_sampled)[names(TBsg17_sampled)=="Freq"] <- "sum"
# join M1 and M2
names(TBsg17)[names(TBsg17)=="sum(Total)"] <- "sum"
TBsg17$Year <- as.character(TBsg17$Year)
TBsg17_final <- bind_rows(TBsg17, TBsg17_sampled)
TBsg17_final <- TBsg17_final %>%
  filter(sum != 0)


# No M2 species
names(TCma17_10)[names(TCma17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TCma17_10_pool <- rep(TCma17_10$Species, TCma17_10$sum)
TCma17_10_pool <- factor(TCma17_10_pool)
set.seed(521)
TCma17_sampled <- as.data.frame(table((sample(TCma17_10_pool, sum(TCma17_10$sum)/5, replace = FALSE, prob = rep(1/sum(TCma17_10$sum), length(TCma17_10_pool))))))
TCma17_sampled$Site.Name <- "Twin Cays Mangrove"
TCma17_sampled$Habitat <- "Mangrove"  
TCma17_sampled$Year <- "2017"
names(TCma17_sampled)[names(TCma17_sampled)=="Var1"] <- "Species"
names(TCma17_sampled)[names(TCma17_sampled)=="Freq"] <- "sum"
# filter zero draws
TCma17_final <- TCma17_sampled %>%
  filter(sum != 0)


# No M2 species
names(TCsg17_10)[names(TCsg17_10)=="sum(Total)"] <- "sum"
# create pool of M1 only fishes to draw from
TCsg17_10_pool <- rep(TCsg17_10$Species, TCsg17_10$sum)
TCsg17_10_pool <- factor(TCsg17_10_pool)
set.seed(682)
TCsg17_sampled <- as.data.frame(table((sample(TCsg17_10_pool, sum(TCsg17_10$sum)/5, replace = FALSE, prob = rep(1/sum(TCsg17_10$sum), length(TCsg17_10_pool))))))
TCsg17_sampled$Site.Name <- "Twin Cays Seagrass"
TCsg17_sampled$Habitat <- "Seagrass"  
TCsg17_sampled$Year <- "2017"
names(TCsg17_sampled)[names(TCsg17_sampled)=="Var1"] <- "Species"
names(TCsg17_sampled)[names(TCsg17_sampled)=="Freq"] <- "sum"
# filter zero draws
TCsg17_final <- TCsg17_sampled %>%
  filter(sum != 0)


######## BIND ALL SITES INTO SINGLE DATASET

Subsamp_RLS <- bind_rows(BGma15_final, BGma16_final, BGma17_final, CBfr15_final, CBfr16_final, CBfr17_final, CBpa15_final, CBpa16_final, CBpa17_final, CBsa15_final, CBsa16_final, CBsa17_final, CBsg15_final, CBsg16_final, CBsg17_final, CUpr15_final, CUpr16_final, CUpr17_final, CUsa15_final, CUsa16_final, CUsa17_final, CUsg15_final, CUsg16_final, CUsg17_final, HRpr15_final, HRpr16_final, HRpr17_final, SRfr15_final, SRfr16_final, SRfr17_final, TBfr15_final, TBfr16_final, TBfr17_final, TBma15_final, TBma16_final, TBma17_final, TBsa15_final, TBsa16_final, TBsa17_final, TBsg17_final, TCma15_final, TCma16_final, TCma17_final, TCsg15_final, TCsg16_final, TCsg17_final)

# write.csv(Subsamp_RLS, "20171120_CBC_RLS_subsampled.csv")


###################################################################################
# END OF SCRIPT                                                                   #
###################################################################################

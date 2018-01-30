###################################################################################
#                                                                                ##
# Biomass Calculations for Carrie Bow RLS Data                                   ##
# Data are current as of 2018-01-30                                              ##
# Data source: Smithsonian Institution                                           ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-01-30                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:
# This script calculates biomass of fish observed on Reef Life Surveys conducted
# at Carrie Bow Cay from 2015-2017. Only considers method 1 (pelagic fishes).

# Required Files (check that script is loading latest version):
# CBCRLS_CORE.csv
# 20180116_RLS_biomass_coefs.csv

# Associated Scripts:
# CBC_Manuscript_Figures.R 

# TO DO

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# SUMMARY STATS                                                                   #   
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2018-01-30  Script Created. 

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse) # plots #data manipulation

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# set working directory
setwd("~/Dropbox (Personal)/Datasets/CBC_MS_2018")

# Import full survey for method 1 biomass calcs
RLS_1 <- read.csv("CBCRLS_CORE.csv")

# Load coefficients
coef <- read.csv("20180116_RLS_biomass_coefs.csv")

###################################################################################
# BIOMASS CALCULATIONS                                                            #
###################################################################################

# Filter only Method 1
RLS_1 <- RLS_1 %>%
  filter(Method == "1")

# Clean up data

# Create unique ID for each survey
RLS_1 <- RLS_1 %>%
  unite_("event", c("Site.Name", "Date"), sep = "_", remove = FALSE)

# make unique ID a factor
RLS_1$event <- as.factor(RLS_1$event)

#remove method 1 with no given abundance
RLS_1 <- RLS_1 %>%
  subset(Total !="0")

# fix incorrect taxa
levels(RLS_1$Species)[levels(RLS_1$Species)=="Sphoeroides spenglerii"] <- "Sphoeroides spengleri"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Xyricthys splendens"] <- "Xyrichtys splendens"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Ctenogonius stigmaturus"] <- "Ctenogobius stigmaturus"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Ophistignathus sp."] <- "Opistognathus sp."
levels(RLS_1$Species)[levels(RLS_1$Species)=="Ctenogonius stigmaturus"] <- "Ctenogobius stigmaturus"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Hyppolytidae"] <- "Hippolytidae"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Paradiplogrammus bairdi"] <- "Callionymus bairdi"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Cryptotoums roseus"] <- "Cryptotomus roseus"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Eucinostoums gula"] <- "Eucinostomus gula"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Atherinidae"] <- "Atherinid spp."
levels(RLS_1$Species)[levels(RLS_1$Species)=="Clupeoid spp."] <- "Clupeidae"
levels(RLS_1$Species)[levels(RLS_1$Species)=="Antherinidae sp."] <- "Atherinid spp."


# remove unknown and freshwater species
RLS_1 <- RLS_1 %>%
  filter(Species != c("Gambusia sp.")) %>%
  filter(Species != c("Scarine spp."))

# remove incidental inverts
RLS_1 <- RLS_1 %>%
  filter(Species != "Hippolytidae") %>%
  filter(Species != "Mithracid spp.")%>%
  filter(Species != "Strombus gigas") %>%
  filter(Species != "Sepioteuthis sepioidea")

# replace NAs with 0 
RLS_1[is.na(RLS_1)] <- 0



############### BIOMASS COEFFICIENTS

# rename species column to match RLS column
names(coef)[names(coef)=="SPECIES_NAME"] <- "Species"

# remove unused columns
coef_small <- coef %>%
  select(Species, A, B, CF)

# match characters to RLS data
RLS_biom <- full_join(RLS_1, coef_small, by="Species")

# identify species not in character list all M1M2
unique(RLS_biom$Species[!(RLS_biom$Species %in% coef$Species)])

write.csv(RLS_biom, "CBC_RLS_fish-biomass_full.csv")

### CALCULATIONS FOR BIOMASS MADE IN EXCEL AND RELOADED INTO R as 
### "20180124_CBC_biom.csv"


#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#
###################################################################################
#                                                                                ##
# Carrie Bow Cay MarineGEO Manuscript Figures                                    ##
# Data are current as of 20171122                                                ##
# Data source: MarineGEO - Tennenbaum Marine Observatories Network - Smithsonian ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 20171122                                                          ##
#                                                                                ##
###################################################################################

# SUMMARY:
# Simplified script derived from RLS_subsampled_CBC_analyses_v2.R, 
# CBC_Squidpops_v2.R, and 20170228_weedpops_CBC.R for creation of figures in       
# MarineGEO Carrie Bow Cay Manuscript.

# Required Files (check that script is loading latest version):
# 20171120_CBC_RLS_subsampled.csv
# Traits_all-species_edit.csv
# 20171122_squidpops_CBC.csv
# 20171120_weedpops.csv
# COPY_20171206_CBC-Benthic-photo-annotations.csv
# 20170209_CBC_Relief.csv

# Associated Scripts:
# RLS_M1M2_v1.R 

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# RICHNESS & ABUNDANCE                                                            #
# TAXONOMIC DIVERSITY                                                             #
# COMMUNITY COMPOSITION                                                           #
# TIME SERIES                                                                     #
# PAIRS DATA                                                                      #
# BENTHIC SPACE                                                                   #
# FACETED FIGURES                                                                 #
#                                                                                 #
###################################################################################

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(ggplot2) # plots
library(plyr) #data manipulation
library(dplyr) # data manipulation
library(tidyr) # data manipulation
library(ggpubr) # multi-plot visualizations
library(vegan) # diversity and MDS plots
library(splitstackshape) # data manipulation
library(viridis) # color palette
library(psych) # pairs analysis

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# Change working directory
# setwd("~/Documents/Git/Carrie Bow Analyses")

############### RLS

# Import RLS survey data
RLS_subsample <- read.csv("20171120_CBC_RLS_subsampled.csv")
RLS_subsample$Year <- as.factor(RLS_subsample$Year)
RLS_subsample$Habitat <- factor(RLS_subsample$Habitat, levels = c("Fore Reef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))
glimpse(RLS_subsample)

############### FISH CHARACTERS

# import character list
characters <- read.csv("Traits_all-species_edit.csv")

# fix level duplication
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="planktivore"] <- "Planktivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="higher carnivore"] <- "Higher carnivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="omnivore"] <- "Omnivore"
levels(characters$Trophic.group)[levels(characters$Trophic.group)=="suspension feeders"] <- "Suspension feeders"
levels(characters$Trophic.group)

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

############### BENTHIC COMPOSITION

benthic <- read.csv("COPY_20171206_CBC-Benthic-photo-annotations.csv")

# glimpse data
glimpse(benthic)

# change year to factor
benthic$Year <- as.factor(benthic$Year)

# remove extraneous columns
short_benthic <- benthic %>%
  select(Name, Site.Name, Year, Habitat, Label)

# Order factor Site.Name by reef type and from North to South
short_benthic$Site.Name <- factor(short_benthic$Site.Name, levels = c("Tobacco Reef", "Carrie Bow Reef", "South Reef", "Carrie Bow Patch Reef", "House Reef", "Curlew Patch Reef"))

# add up occurences 
short_benthic$amount <- rep(1)
short_benthic <- short_benthic %>%
  group_by(Name, Site.Name, Year, Habitat, Label, amount) %>%
  summarise(sum(amount))

# change column name
names(short_benthic)[names(short_benthic)=="sum(amount)"] <- "sum"

# spread data
short_benthic <- short_benthic %>%
  spread(Label, sum) 

# replace NAs with 0 
short_benthic[is.na(short_benthic)] <- 0

# remove "amount" column, and null data points/unknowns  
short_benthic <- short_benthic %>%
  select(c(-TAPE, -amount, -Unk, -SHAD))

# create biotic-only dataset for MDS anlayses
biotic_benthic <- short_benthic %>%
  select(-LIME, -Sand, -LSUB_RUB)

# create biotic-only dataset for diversity anlayses
benthic_rich <- short_benthic 

# create dataset for space calcs
space_benthic <- short_benthic

############### REEF RELIEF

# read in relief data
relief <- read.csv("20170209_CBC_Relief.csv")

# change year to factor
relief$Year <- as.factor(relief$Year)

# Order factor Site.Name by reef type and from North to South
relief$Site.Name <- factor(relief$Site.Name, levels = c("Tobacco Reef", "Carrie Bow Reef", "South Reef", "Carrie Bow Patch Reef", "House Reef", "Curlew Patch Reef"))

# remove unused columns
relief <- relief %>%
  select(c(Habitat, Site.Name, Year, Relief..m.))

# create column for joining to other data
relief <- relief %>%
  unite_("Site.Time", c("Site.Name", "Year"), sep = "_", remove = FALSE)


###################################################################################
# RICHNESS & ABUNDANCE                                                            #
###################################################################################

############### FISH RICHNESS

#set dataset to analyse
#year, method
RLS_spread <- RLS_subsample %>%
  spread(Species, sum, fill = 0)

RLS_rich <- RLS_spread %>%
  group_by(Site.Name, Habitat, Year) %>%
  summarise_all(funs(sum))

RLS_rich <- RLS_rich %>%
  select(-X)

pool <- RLS_rich[,4:123]
sitetime <- RLS_rich[,c(1:3)]
sitetime <- sitetime %>%
  unite_("Site.Time", c("Site.Name", "Year"), sep = "_", remove = FALSE)
attach(sitetime)
pool <- RLS_rich[,4:123] %>%
  specpool(Site.Time)

############### TOTAL FISH ABUNDANCE

detach(package:plyr)
Abundance <- RLS_subsample %>%
  group_by(Site.Name, Year) %>%
  summarise(log10(sum(sum)))

Abundance$Habitat <- c("Mangrove", "Mangrove", "Mangrove", "Patch Reef", "Patch Reef", "Patch Reef", "Fore Reef", "Fore Reef",  "Fore Reef", "Sand", "Sand", "Sand", "Seagrass", "Seagrass", "Seagrass", "Patch Reef", "Patch Reef", "Patch Reef", "Sand",  "Sand", "Sand", "Seagrass", "Seagrass", "Seagrass", "Patch Reef", "Patch Reef", "Patch Reef", "Fore Reef", "Fore Reef", "Fore Reef", "Mangrove",  "Mangrove", "Mangrove", "Fore Reef", "Fore Reef", "Fore Reef", "Sand", "Sand", "Sand", "Seagrass",  "Mangrove", "Mangrove", "Mangrove", "Seagrass", "Seagrass", "Seagrass")

Abundance$Habitat <- factor(Abundance$Habitat, levels = c("Fore Reef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

Abundance$Richness <- pool$Species
names(Abundance)[names(Abundance)=="log10(sum(sum))"] <- "Abun"

###################################################################################
# TAXONOMIC DIVERSITY                                                             #
###################################################################################

#year, method
RLS_div <- RLS_subsample %>%
  spread(Species, sum, fill = 0)
sp_simp <- RLS_div[,c(2:124)]
sp_simp <- sp_simp %>%
  group_by(Site.Name, Year, Habitat) %>%
  summarise_all(funs(sum))
y <- sp_simp[,4:123]
z <- sp_simp[,c(1:3)]

attach(z)

############### FISH SPECIES DIVERSITY

# simpson diversity
sp_simpson <- diversity(y, index = "simpson")
z$sp_simpson <- sp_simpson

############### FISH FAMILY DIVERSITY

RLS_func <- RLS_subsample

# match characters to RLS data
names(characters)[names(characters)=="CURRENT_TAXONOMIC_NAME"] <- "Species"
RLS_func <- left_join(RLS_func, characters, by="Species")

# identify species not in character list
unique(RLS_func$Species[!(RLS_func$Species %in% characters$Species)])

# remove extraneous columns for family analysis
RLS_fam <- RLS_func %>%
  select(Site.Name, Habitat.x, Year, Species, sum, Family, Class)

# remove 0 occurrences
RLS_fam <- RLS_fam %>%
  filter(sum != "0")

# sum familial groups for each site 
Fam_anal <- RLS_fam %>%
  group_by(Site.Name, Habitat.x, Year, Family) %>%
  summarise(sum(sum))

# change column name
names(Fam_anal)[names(Fam_anal)=="sum(sum)"] <- "total"

# spread for diversity analysis
Fam_anal <- Fam_anal %>%
  spread(Family, total, fill = 0)

# simpson diversity of familial groups
fam_simp <- Fam_anal[,5:48]
fam_cats <- Fam_anal[,1:3]
names(fam_cats)[names(fam_cats)=="Habitat.x"] <- "Habitat"
attach(fam_cats)

simpson <- diversity(fam_simp, index = "simpson")
fam_cats$simpson <- simpson

############### HERBIVORE ABUNDANCE & RICHNESS

Herbs <- RLS_func %>%
  select(Site.Name, Year, Habitat.x, Species, sum, Trophic.group) %>%
  filter(Trophic.group %in% c("Algal farmer", "Browsing herbivore", 
                              "Scraping herbivore"))
  
Herb.abun <- Herbs %>%
  group_by(Site.Name, Year, Habitat.x) %>%
  summarise(sum(sum))
names(Herb.abun)[names(Herb.abun)=="sum(sum)"] <- "total"

Herb.spread <- Herbs %>%
  group_by(Site.Name, Year, Habitat.x, Species) %>%
  summarise(sum(sum))
names(Herb.spread)[names(Herb.spread)=="sum(sum)"] <- "total"

Herb.spread <- Herb.spread %>%
  spread(Species, total, fill = 0)

herb.pool <- Herb.spread[,4:27]
herb.sitetime <- Herb.spread[,c(1:3)]
herb.sitetime <- herb.sitetime %>%
  unite_("Site.Time", c("Site.Name", "Year"), sep = "_", remove = FALSE)
attach(herb.sitetime)
herb.pool <- Herb.spread[,4:27] %>%
  specpool(Site.Time)

############### BENTHIC DIVERSITY

# add simpson diversity values in column
benthic_rich$simpson <- diversity(benthic_rich[,5:55], index = "simpson")

###################################################################################
# COMMUNITY COMPOSITION                                                           #
###################################################################################

############### FAMILIAL COMPOSITION

# remove extraneous columns - use this dataframe for family stacked plots
RLS_clean <- RLS_func %>%
  select(Site.Name, Habitat.x, Trophic.group, Family, sum, Year)
#RLS_clean <- RLS_clean[!(RLS_clean$Trophic.group == "NA"),]

hab_fam <- RLS_clean %>%
  group_by(Habitat.x, Family) %>%
  summarise(sum(sum))
# rename column
names(hab_fam)[names(hab_fam)=="sum(sum)"] <- "total"

### Calculate species for each hab that are < 5% total abundance
#fr <- hab_fam %>%
#  filter(Habitat.x == "Fore Reef") %>%
#  summarise(sum(total))
# 2431* 0.05 = 121.55 5%
# 2431* 0.03 = 72.93 3%
# 2431* 0.01 = 24.31 1%
#pr <- hab_fam %>%
#  filter(Habitat.x == "Patch Reef") %>%
#  summarise(sum(total))
# 2870* 0.05 = 143.5 5%
# 2870* 0.03 = 86.1 3%
# 2870* 0.01 = 28.7 1%
#ma <- hab_fam %>%
#  filter(Habitat.x == "Mangrove") %>%
#  summarise(sum(total))
# 6579* 0.05 = 328.95 5%
# 6579* 0.03 = 197.37 3%
# 6579* 0.01 = 65.79 1%
#sa <- hab_fam %>%
#  filter(Habitat.x == "Sand") %>%
#  summarise(sum(total))
# 569* 0.05 = 28.45 5%
# 569* 0.03 = 17.07 3%
# 569* 0.01 = 5.69 1%
# sg <- hab_fam %>%
#  filter(Habitat.x == "Seagrass") %>%
#  summarise(sum(total))
# 1369* 0.05 = 68.45 5%
# 1369* 0.03 = 41.07 3%
# 1369* 0.01 = 13.69 1%

# Lump families <5% of total abundance per habitat as "Other"
habfamlev <- levels(hab_fam$Family)
levels(hab_fam$Family) = c(habfamlev, "Other")
hab_fam[hab_fam$Habitat.x == "Fore Reef" & hab_fam$total <122, 2] <- "Other"
hab_fam[hab_fam$Habitat.x == "Patch Reef" & hab_fam$total <144, 2] <- "Other"
hab_fam[hab_fam$Habitat.x == "Mangrove" & hab_fam$total <329, 2] <- "Other"
hab_fam[hab_fam$Habitat.x == "Sand" & hab_fam$total <29, 2] <- "Other"
hab_fam[hab_fam$Habitat.x == "Seagrass" & hab_fam$total <69, 2] <- "Other"


freq_fam <- hab_fam %>% 
  expandRows("total")
# reorder factor levels
freq_fam$Habitat.x <- factor(freq_fam$Habitat.x, levels = c("Fore Reef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))

totals <- RLS_clean %>%
  group_by(Habitat.x) %>%
  summarize(total = sum(sum))

############### TROPHIC COMPOSITION

# sum tropic groups per habitat type
hab_troph <- RLS_clean %>%
  group_by(Habitat.x, Trophic.group) %>%
  summarise(sum(sum))
# rename column
names(hab_troph)[names(hab_troph)=="sum(sum)"] <- "total"

freq_troph <- hab_troph %>%
  expandRows("total")
# reorder factor levels
freq_troph$Habitat.x <- factor(freq_troph$Habitat.x, levels = c("Fore Reef", "Patch Reef", "Mangrove", "Seagrass", "Sand"))
# remove rows with NA
freq_troph <- na.omit(freq_troph)
# remove blank rows
freq_troph <-freq_troph[!(freq_troph$Trophic.group ==""),]

############### TAXONOMIC MDS

# spread RLS data into community matrix
RLS_full_mat <- RLS_subsample %>%
  spread(Species, sum, fill = 0)
RLS_mat <- RLS_full_mat[,c(2:140)]
RLS_mat <- RLS_mat %>%
  group_by(Site.Name, Habitat, Year) %>%
  summarise_all(funs(sum))

# run matrix 
RLS_mat_run <- RLS_mat[,4:139]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
plot_data_tax <- data.frame(RLS_mds_points, RLS_mat[,1:3])
library(plyr)
chulls_tax <- ddply(plot_data_tax, .(Habitat), function(df) df[chull(df$MDS1, df$MDS2), ])
detach(package:plyr)

############### TROPHIC MDS

# spread RLS data into community matrix
RLS_troph_mat <- RLS_func %>%
  spread(Trophic.group, sum, fill = 0)
str(RLS_troph_mat, list.len=ncol(RLS_troph_mat))
RLS_mat <- RLS_troph_mat[,c(2:4, 17:25)]
names(RLS_mat)[names(RLS_mat)=="Habitat.x"] <- "Habitat"
RLS_mat <- RLS_mat %>%
  group_by(Site.Name, Habitat, Year) %>%
  summarise_all(funs(sum))

# run matrix 
RLS_mat_run <- RLS_mat[,4:12]

# MDS of community
RLS_mds <- metaMDS(RLS_mat_run)
RLS_mds_points <- RLS_mds$points
RLS_mds_points <- data.frame(RLS_mds_points)
plot_data_tro <- data.frame(RLS_mds_points, RLS_mat[,1:3])
library(plyr)
chulls_tro <- ddply(plot_data_tro, .(Habitat), function(df) df[chull(df$MDS1, df$MDS2), ])
detach(package:plyr)

############### BENTHIC MDS

# run matrix 
benthic_mat_run <- biotic_benthic[,5:52]

# MDS of benthic community
ben_mds <- metaMDS(benthic_mat_run)
ben_mds_points <- ben_mds$points
ben_mds_points <- data.frame(ben_mds_points)
plot_data_ben <- data.frame(ben_mds_points, biotic_benthic[,1:4])

# hulls by habitat
library(plyr)
#chulls_ben <- ddply(plot_data_ben, .(Habitat), function(df) df[chull(df$MDS1, df$MDS2), ])

# hulls by year
chulls_ben_year <- ddply(plot_data_ben, .(Year), function(df) df[chull(df$MDS1, df$MDS2), ])

# hulls by site
chulls_ben_site <- ddply(plot_data_ben, .(Site.Name), function(df) df[chull(df$MDS1, df$MDS2), ])
detach(package:plyr)

###################################################################################
# TIME SERIES                                                                     #
###################################################################################

# sum total macrophtye lost
Algae <- weedpops %>%
  group_by(year, Line_ID, habitat) %>%
  summarise(sum(detachment.24hr))
names(Algae)[names(Algae)=="sum(detachment.24hr)"] <- "detachment"

###################################################################################
# PAIRS DATA                                                                      #
###################################################################################

# bind herbivore data sets together
Total.herbs <- bind_cols(herb.pool, herb.sitetime)
Herb.abun <- Herb.abun %>%
  unite_("Site.Time", c("Site.Name", "Year"), sep = "_", remove = FALSE)
Total.herbs <- left_join(Total.herbs, Herb.abun, by="Site.Time")

# create common index for community richness/abundance data to join to herbs
Abundance <- Abundance %>%
  unite_("Site.Time", c("Site.Name", "Year"), sep = "_", remove = FALSE)

# create common index for squidpop data to join with others
squidpops <- squidpops %>%
  unite_("Site.Time", c("Site", "Year"), sep = "_", remove = FALSE)

# create common index for weedpop data to join with others
weedpops <- weedpops %>%
  unite_("Site.Time", c("site", "year"), sep = "_", remove = FALSE)
weedpop.wide <- weedpops %>%
  select(Line_ID, Site.Time, habitat, item, detachment.24hr)
weedpop.wide <- weedpop.wide %>%
  spread(item, detachment.24hr, fill = 0)
weedpop.wide$Tot.algae <- rowSums(weedpop.wide[,4:8])

# join abundance all with herb and squidpop data
All_pairs <- left_join(Abundance, Total.herbs, by="Site.Time")
All_pairs <- left_join(All_pairs, squidpops, by="Site.Time")
All_pairs <- left_join(All_pairs, weedpop.wide, by="Site.Time")
names(All_pairs)[names(All_pairs)=="Species"] <- "Herb.rich"
names(All_pairs)[names(All_pairs)=="total"] <- "Herb.abun"
names(All_pairs)[names(All_pairs)=="Abun"] <- "Tot.log.abun"
names(All_pairs)[names(All_pairs)=="Richness"] <- "Tot.rich"
names(All_pairs)[names(All_pairs)=="Detachment.1.hour"] <- "Squidpops"
All_pairs <- All_pairs %>%
  select(Site.Time, Habitat, Tot.log.abun, Tot.rich, Herb.abun, Herb.rich, 
         Squidpops, Tot.algae, acanthophora, dictyota, halimeda, sargassum, 
         thalassia)

# calculate mean rugosity per site/time
relief_mean <- relief %>%
  group_by(Site.Time) %>%
  summarize(mean(Relief..m.))

# join relief data with all other data
All_pairs <- left_join(All_pairs, relief_mean, by="Site.Time")

# rename calculated column
names(All_pairs)[names(All_pairs)=="mean(Relief..m.)"] <- "mean.relief"

###################################################################################
# BENTHIC SPACE                                                                   #
###################################################################################

# add total number of points per photo
space_benthic$points <- rowSums(space_benthic[,5:55])

# add total number of biotic points per photo
space_benthic$empty <- rowSums(space_benthic[,c(8, 22, 24, 43)])

# calc occupied space
space_benthic$occupied <- space_benthic$points - space_benthic$empty
space_benthic$empty_prop <- space_benthic$empty/space_benthic$points
empty_space <- space_benthic %>%
  group_by(Name, Site.Name, Year, Habitat) %>%
  summarise(median(empty_prop))

# add total number of macroalgal points per photo
space_benthic$macroalgae <- rowSums(space_benthic[,c(12,13,21,23,25,44,54,55)])
space_benthic$algal_prop <- space_benthic$macroalgae/space_benthic$points
algal_cover <- space_benthic %>%
  group_by(Name, Site.Name, Year, Habitat, algal_prop) %>%
  summarise()

#rename column
names(algal_cover)[names(algal_cover)=="algal_prop"] <- "total"

###################################################################################
# FACETED FIGURES                                                                 #
###################################################################################

habcols <- c("#D55E00", "#CC79A7", "#e79f00", "#009E73", "#9ad0f3")
habcols_long <- as.character(c("Fore Reef"="#D55E00", "Patch Reef"="#CC79A7", "Mangrove"="#e79f00", "Seagrass"="#009E73", "Sand"="#9ad0f3"))
habScale <- scale_colour_manual(name="Habitats", values=habcols_long)

############### Figure 1

# 1A Fish Abundance
abun.plot <- ggplot(Abundance, aes(x= Habitat, y = Abun, fill = Habitat)) + 
  geom_boxplot() + 
  scale_fill_manual(values = habcols) +
  geom_point(aes(x=Habitat, y=Abun, shape=Habitat)) + 
  scale_y_continuous(limits = c(0,4)) + 
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Log Abundance")
#abun.plot

# 1B Fish Species Richness
rich.plot <- ggplot(Abundance, aes(x= Habitat, y = Richness, fill = Habitat)) + 
  geom_boxplot()  + 
  scale_fill_manual(values = habcols) +
  geom_point(aes(x=Habitat, y=Richness, shape=Habitat)) + 
  scale_y_continuous(limits = c(0,35)) + 
  theme_minimal() + 
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Number of Species")
#rich.plot

# 1C Fish Species Simpson Diversity
spdiv.box <- ggplot(sp_simp, aes(x=Habitat, y=sp_simpson, fill = Habitat)) +
  geom_boxplot() +
  scale_fill_manual(values = habcols) +
  geom_point(aes(x=Habitat, y=sp_simpson, shape = Habitat)) + 
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Simpson Value (species)")
#spdiv.box

# 1D Fish Familial Simpson Diversity
fadiv.box <- ggplot(fam_simp, aes(x=Habitat, y=simpson, fill = Habitat)) + 
  geom_boxplot() + 
  scale_fill_manual(values = habcols) + 
  geom_point(aes(x=Habitat, y=simpson, shape = Habitat)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Simpson Value (family)")
#fadiv.box

Figure1 <- ggarrange(abun.plot, rich.plot, spdiv.box, fadiv.box, 
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2,
                     common.legend = TRUE, legend = "right")
annotate_figure(Figure1, bottom = text_grob("Figure 1: Observed measures of A) fish abundance, B) fish species richness, C) fish species \n simpson diversity, and D) fish family simpson diversity across five habitat types sampled \n annually from 2015-17 around Carrie Bow Cay, Belize.", size = 10))

# best size: ~630x700

############### FIGURE 2

# 2A Familial Composition
fam_plot <- ggplot(freq_fam, aes(Habitat.x)) + 
  geom_bar(aes(fill = Family), position = "fill") + 
  theme_minimal() + 
  scale_fill_viridis(discrete=TRUE) +
  labs(x="", y="Proportion") + 
  geom_text(data = totals, aes(x = c(1:5), y = 1.05, label = total)) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)) 
#fam_plot

# 2B Trophic Group Composition
troph_plot <- ggplot(freq_troph, aes(Habitat.x), na.rm = TRUE) + 
  geom_bar(aes(fill = Trophic.group), position = "fill") + 
  theme_minimal() + 
  scale_fill_viridis(discrete = TRUE, option = "C", name ="Trophic Group") +
  labs(x="", y="Proportion") + 
  geom_text(data = totals, aes(x = c(1:5), y = 1.05, label = total)) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)) 
#troph_plot

# 2C Taxonomic nMDS
tax_MDS <- ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, pch = Year, 
                                     color = Habitat)) +  
  theme_minimal() +
  geom_point(size = 4) +
  geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) +
  habScale + 
  rremove("legend")
#tax_MDS

# 2D Trophic nMDS
tro_MDS <- ggplot(plot_data_tro, aes(x=MDS1, y=MDS2, pch = Year, 
                                     color = Habitat)) + 
  theme_minimal() +
  geom_point(size = 4) + 
  geom_polygon(data=chulls_tro, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) +
  habScale
#tro_MDS

Figure2 <- ggarrange(fam_plot, troph_plot, tax_MDS, tro_MDS, 
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2,
                     legend = NULL)
annotate_figure(Figure2, bottom = text_grob("Figure 2: Calcuated A) family composition, B) trophic groups, C) species composition, and D) trophic composition of fish communities \n across five habitat types sampled annually from 2015-17 around Carrie Bow Cay, Belize. Total number of fishes included in \n proportional measures shown above each bar", size = 10))

# best size: ~820x820

############### FIGURE 3

# 3A One Hour Squidpops
one.hour.box <- ggplot(squidpops, aes(x=Habitat, y=Proportion.Missing.1.hour, 
                                      fill=Habitat)) + 
  geom_boxplot()  + 
  scale_fill_manual(values = habcols) +
  geom_point(aes(x=Habitat, y=Proportion.Missing.1.hour, shape=Habitat)) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme_minimal() + 
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Proportion Missing")
#one.hour.box

# 3B Twenty-four Hour Weedpops
tfhour.weedpop <- ggplot(weedpops, aes(x=habitat, y=prop.24hr, fill = item)) +  
  geom_boxplot(position = "dodge") + 
  scale_fill_viridis(discrete = TRUE, option = "A", name ="Macrophyte", 
                     begin = 0.4) +
  geom_line(aes(group = habitat)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)) +
  labs(x="", y="Proportion Missing")
#tfhour.weedpop

Figure3 <- ggarrange(one.hour.box, tfhour.weedpop, 
                     labels = c("A", "B"),
                     ncol = 1, nrow = 2)
annotate_figure(Figure3, bottom = text_grob("Figure 3: The A) squidpops at one hour and B) weedpops at twenty-four hours.", size = 10))

# best size: ~570x620

############### FIGURE 4

# Figure 4A timeseries of richness
richness.time <- ggplot(Abundance, aes(x=Year, y=Richness, group = Habitat)) + 
  geom_point(size = 2, position = position_jitter(w = 0.05, h = 0), 
             aes(color = Habitat)) + 
  stat_summary(fun.y=mean, geom="line", size = 1, aes(group=Habitat, 
                                                      colour = Habitat))  +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Number of Species") +
  habScale
#richness.time

# Figure 4B timeseries of abundance
abun.time <- ggplot(Abundance, aes(x=Year, y=Abun, group = Habitat)) + 
  geom_point(size = 2, position = position_jitter(w = 0.05, h = 0), 
             aes(color = Habitat)) + 
  stat_summary(fun.y=mean, geom="line", size = 1, aes(group=Habitat, 
                                                      colour = Habitat))  +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Log Abundance") +
  habScale
#abun.time

# Figure 4C timeseries of total squidpop consumption
squidpop.time <- ggplot(squidpops, aes(x=Year, y=Detachment.1.hour, 
                                       group = Habitat)) + 
  geom_point(size = 2, position = position_jitter(w = 0.05, h = 0), 
             aes(color = Habitat)) + 
  stat_summary(fun.y=mean, geom="line", size = 1, aes(group=Habitat, 
                                                      colour = Habitat))  +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Bait Loss") +
  habScale
#squidpop.time

# Figure 4D timeseries of total algal consumption
algae.time <- ggplot(Algae, aes(x=year, y=detachment, group = habitat)) + 
  geom_point(size = 2, position = position_jitter(w = 0.05, h = 0), 
             aes(color = habitat)) + 
  stat_summary(fun.y=mean, geom="line", size = 1, aes(group=habitat, 
                                                      colour = habitat))  +
  theme_minimal() +
  labs(x="", y="Algal Loss") +
  habScale
#algae.time

Figure4 <- ggarrange(richness.time, abun.time, squidpop.time, algae.time, 
                     labels = c("A", "B", "C", "D"),
                     ncol = 1, nrow = 4,
                     common.legend = TRUE, legend = "right")
annotate_figure(Figure4, bottom = text_grob("Figure 4: Calcuated A) species richness, B) log abundance, C) total  \n squidpop bait loss, and D) total algal loss across five habitat types \n sampled annually from 2015-17 around Carrie Bow Cay, Belize.", size = 10))

# best size: ~420x900

############### FIGURE 5

# Figure 5A pairs analysis of all habitat types

all_pairplot <- pairs.panels(All_pairs[5:15],bg=habcols[All_pairs$Habitat.x], 
                             scale = TRUE, pch=21, ellipses = FALSE)
# best size: ~1150x850

# Figure 5B pairs for reef only

reef_pairplot <- pairs.panels(All_pairs[c(4:5, 7:8, 19:20, 29:30, 32:33, 38:39),
                                        c(5:16)], bg=habcols[All_pairs$Habitat.x],
                              scale = FALSE, pch=21, ellipses = FALSE)
# best size: ~1150x850

Figure5 <- ggarrange(all_pairplot, reef_pairplot, ncol = 1, nrow = 2)
annotate_figure(Figure5, bottom = text_grob("Figure 5: Pairs plot of stuff."))

############### FIGURE 6

# Figure 6A Benthic MDS Year
benthic_MDS_Year <- ggplot(plot_data_ben, aes(x=MDS1, y=MDS2, pch = Habitat, 
                                         color = Year)) + 
  scale_color_viridis(discrete = TRUE, option = "viridis", begin = 0.3, end = 0.7) +
  theme_minimal() + 
  geom_point(size = 4) + 
  geom_polygon(data=chulls_ben_year, aes(x=MDS1, y=MDS2, group=Year), fill=NA) 
#benthic_MDS_Year

# Figure 6x Benthic MDS Site (sites all overlap)
#benthic_MDS_Site <- ggplot(plot_data_ben, aes(x=MDS1, y=MDS2, pch = Year, 
#                                         color = Site.Name)) + 
#  scale_color_viridis(discrete = TRUE, option = "viridis", begin = 0.1, end = 0.9) +
#  theme_minimal() + 
#  geom_point(size = 4) + 
#  geom_polygon(data=chulls_ben_site, aes(x=MDS1, y=MDS2, group=Site.Name), fill=NA) 
#benthic_MDS_Site

# Figure 6B Unoccupied Space
benspace.box <- ggplot(space_benthic, aes(x = Site.Name, y = 1-(occupied/points), fill = Year)) + 
  geom_boxplot(aes(fill = Year, group = interaction(Site.Name, Year))) + 
  scale_fill_viridis(discrete = TRUE, option = "viridis", begin = 0.3, end = 0.7) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Proportion Space")
#benspace.box

# Figure 6C Macroalgal Abundance
benalgae.box <- ggplot(algal_cover, aes(x = Site.Name, y = total, fill = Year)) + 
  geom_boxplot(aes(fill = Year, group = interaction(Site.Name, Year))) + 
  scale_fill_viridis(discrete = TRUE, option = "viridis", begin = 0.3, end = 0.7) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Proportion Algae")
#benalgae.box

# Figure 6D Benthic Simpson Diversity
bendiv.box <- ggplot(benthic_rich, aes(x=Site.Name, y=simpson, fill = Year)) + 
  geom_boxplot(aes(fill = Year, group = interaction(Site.Name, Year))) + 
  scale_fill_viridis(discrete = TRUE, option = "viridis", begin = 0.3, end = 0.7) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Simpson Value") 
#bendiv.box

Figure6 <- ggarrange(benthic_MDS_Year, ggarrange(benspace.box, benalgae.box, bendiv.box, 
                     labels = c("B", "C", "D"),
                     nrow = 3, legend = "none"), labels = "A", ncol = 2,
                     common.legend = TRUE, legend = "right")
annotate_figure(Figure6, bottom = text_grob("Figure 6: Benthic reef community A) community composition, B) unoccupied space, C) algal abundance, and \n D) Simpson diversity across six fore- and patch reefs sampled in 2015 and 2017 around Carrie Bow Cay, Belize.", size = 10))

# best size: ~900x630

############### FIGURE 7

# Figure 7 Benthic Rugosity of Reefs

rug.box <- ggplot(relief, aes(x=Site.Name, y=Relief..m., fill = Year)) + 
  geom_boxplot(aes(fill = Year, group = interaction(Site.Name, Year))) + 
  scale_fill_viridis(discrete = TRUE, option = "viridis", begin = 0.3, end = 0.7) +
  theme_minimal() +
  labs(x="", y="Rugosity (m)", caption = "Figure 7: Mean rugosity of reef sites around Carrie Bow Cay, Belize by year.")
rug.box

# best size: ~800x550

############### FIGURE 8

# Figure 8A squidpops and fish abundance
sp_abun <- ggplot(All_pairs, aes(y=Squidpops, x=Tot.log.abun)) + 
  geom_point(size = 2, aes(color = Habitat.x)) +
  geom_ribbon(stat = 'smooth', method = "lm", se = TRUE, alpha=0.05, aes(color = NULL)) +
  geom_line(stat = 'smooth', method = lm) +
  theme_minimal() +
  labs(y="Total squidpops consumed", x="Log abundance of fish") +
  habScale
#sp_abun

# Figure 8B squidpops and fish richness
sp_rich <- ggplot(All_pairs, aes(y=Squidpops, x=Tot.rich)) + 
  geom_point(size = 2, aes(color = Habitat.x)) +
  geom_ribbon(stat = 'smooth', method = "lm", se = TRUE, alpha=0.05, aes(color = NULL)) +
  geom_line(stat = 'smooth', method = lm) +
  theme_minimal() +
  labs(y="Total squidpops consumed", x="Fish species richness") +
  habScale
#sp_rich

# Figure 8C algae and fish abundance
alg_abun <- ggplot(All_pairs, aes(y=Tot.algae, x=Tot.log.abun)) + 
  geom_point(size = 2, aes(color = Habitat.x)) +
  geom_ribbon(stat = 'smooth', method = "lm", se = TRUE, alpha=0.05, aes(color = NULL)) +
  geom_line(stat = 'smooth', method = lm) +
  theme_minimal() +
  labs(y="Total algae consumed", x="Log abundance of fish") +
  habScale
#alg_abun


# Figure 8D algae and fish richness
alg_rich <- ggplot(All_pairs, aes(y=Tot.algae, x=Tot.rich)) + 
  geom_point(size = 2, aes(color = Habitat.x)) +
  geom_ribbon(stat = 'smooth', method = "lm", se = TRUE, alpha=0.05, aes(color = NULL)) +
  geom_line(stat = 'smooth', method = lm) +
  theme_minimal() +
  labs(y="Total algae consumed", x="Fish species richness") +
  habScale
#alg_rich

Figure8 <- ggarrange(sp_abun, sp_rich, alg_abun, alg_rich, 
                     labels = c("A","B","C","D"), 
                     ncol = 2,
                     nrow = 2,
                     common.legend = TRUE, legend = "right")
annotate_figure(Figure8, bottom = text_grob("Figure 8: Relationship between squidpops consumed and A) log abundance of fish and B) fish species richness, and total macrophytes \n consumed and C) log abundance of fish and D) fish species richness across 5 habitat types and 3 years around Carrie Bow Cay, Belize.", size = 10))

# best size: ~850x700

###################################################################################
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
###################################################################################



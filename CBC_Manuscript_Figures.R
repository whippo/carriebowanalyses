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

# Associated Scripts:
# RLS_M1M2_v1.R 

# TO DO

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
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
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# Script created from 20171121_CBC_Manuscript_Figures.R and migrated to Git
# 20171206 Added benthic photo MDS analysis
# 20171121 Figures 3-4 created, habitats reordered, saved from 
#          201171117_CBC_Manuscript_Figures.R
# 20171120 Figure 2 created, new colors added, updated data used
# 20171117 Script created from RLS_subsampled_CBC_analyses_v2.R by Ross Whippo

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
setwd("~/Dropbox (Personal)/R_Scripts/CBC 2017/Figures")

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
  select(Site.Name, Year, Habitat, Label)

# add up occurences 
short_benthic$amount <- rep(1)
short_benthic <- short_benthic %>%
  group_by(Site.Name, Year, Habitat, Label, amount) %>%
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

###################################################################################
# RICHNESS & ABUNDANCE                                                            #
###################################################################################

############### RICHNESS

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

############### TOTAL ABUNDANCE

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
benthic_rich$simpson <- diversity(benthic_rich[,4:54], index = "simpson")

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
benthic_mat_run <- biotic_benthic[,4:51]

# MDS of benthic community
ben_mds <- metaMDS(benthic_mat_run)
ben_mds_points <- ben_mds$points
ben_mds_points <- data.frame(ben_mds_points)
plot_data_ben <- data.frame(ben_mds_points, biotic_benthic[,1:3])
library(plyr)
chulls_ben <- ddply(plot_data_ben, .(Habitat), function(df) df[chull(df$MDS1, df$MDS2), ])
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

###################################################################################
# BENTHIC SPACE                                                                   #
###################################################################################

# add total number of points per photo
space_benthic$points <- rowSums(space_benthic[,4:55])

# add total number of biotic points per photo
space_benthic$occupied <- rowSums(space_benthic[,c(4:6, 8:20, 22, 24:41, 43:55)])

# calc empty space
space_benthic$empty <- space_benthic$points - space_benthic$occupied
space_benthic$empty_prop <- space_benthic$empty/space_benthic$points
empty_space <- space_benthic %>%
  group_by(site, date) %>%
  summarise(median(empty_prop))

# add total number of macroalgal points per photo
space_benthic$macroalgae <- rowSums(space_benthic[,c(11,12,20,22,24,43,53,54)])
space_benthic$algal_prop <- space_benthic$macroalgae/space_benthic$points
algal_cover <- space_benthic %>%
  group_by(site,date) %>%
  summarise(median(algal_prop))

space_benthic$site <- factor(space_benthic$site, levels = c("CBC-Lagoon-Reef", "CBC-House-Reef", "CBC-Curlew-Patch-Reef", "CBC-Tobacco-Reef", "CBC-Central-Reef", "CBC-South-Reef"))

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

all_pairplot <- pairs.panels(All_pairs[5:15],bg=habcols[All_pairs$Habitat.x], scale = TRUE,
             pch=21, ellipses = FALSE)

Figure5 <- ggarrange(all_pairplot,
                     ncol = 1, nrow = 1)
annotate_figure(Figure5, bottom = text_grob("Figure 5: Pairs plot of stuff."))

############### FIGURE 6

# Figure 6A Benthic MDS
benthic_MDS <- ggplot(plot_data_ben, aes(x=MDS1, y=MDS2, pch = Year, 
                                     color = Habitat)) + 
  theme_minimal() +
  geom_point(size = 4) + 
  geom_polygon(data=chulls_ben, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) +
  habScale
benthic_MDS

# Figure 6B Benthic Simpson Diversity
ggplot(benthic_rich, aes(site, simpson)) + geom_boxplot() +
  labs(title="Benthic Simpson Diverisy", x="site", y="simpson")

# 1D Fish Familial Simpson Diversity
bendiv.box <- ggplot(benthic_rich, aes(x=Habitat, y=simpson, fill = Habitat)) + 
  geom_boxplot() + 
  scale_fill_manual(values = habcols) + 
  geom_point(aes(x=Habitat, y=simpson, shape = Habitat)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(x="", y="Simpson Value")
bendiv.box


#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#




# This script gives description of the dataset


## Functions--------------------------------------------
library(metafor)
library(dplyr)
library(ggplot2)
library(tidyr)

countstudies <- 
  function(dat, ...){
    return(dat %>% 
             group_by(.dots = lazyeval::lazy_dots(...)) %>% 
             summarize(no.stu = countlevels(ID), no.obs = n()))
  }

'%notin%' <- function(x,y)!('%in%'(x,y))

countlevels <- function(x){
  return(length(levels(as.factor(x))))
}

## Data --------------------------------------------
dat <- read.csv("Output/EffectSizes.csv")
dat$ObsID <- factor(1 : nrow(dat))

### remove NAs-------------------------
dat <- dat[!is.na(dat$ES2),]
dat <- dat[!is.na(dat$VarES2),]

## Description ------------------------------------

# How many studies? 
length(unique(dat$ID))

# How many cases? 
nrow(dat)

# How many pesticides? 
ActiveIngredients <- dat %>% 
  mutate(ActiveIngredients = strsplit(as.character(PollutionNameCorrected), "; ")) %>%
  unnest(ActiveIngredients) %>% select(ActiveIngredients)

nrow(unique(ActiveIngredients)) # 86 substances in total

table(ActiveIngredients)

# Which pesticide categories are reprez?
countstudies(dat, PesticideType)

# what are pesticides most covered? 
PestMost <- data.frame(countstudies(dat, PesticideType, PollutionNameCorrected))
PestMost_stud <- data.frame(PestMost)[order(PestMost$no.stu, decreasing = TRUE),]
head(PestMost_stud, n = 15)

PestMost_obs <- data.frame(PestMost)[order(PestMost$no.obs, decreasing = TRUE),]
head(PestMost_obs, n = 15)

hist(PestMost$no.stu)

# what are pesticides most represented (including in mixtures)
test <- dat %>% mutate(ActiveIngredients = strsplit(as.character(PollutionNameCorrected), "; ")) %>%
  unnest(ActiveIngredients)

PestMost2 <- data.frame(countstudies(test,  ActiveIngredients))
PestMost_stud <- data.frame(PestMost2)[order(PestMost2$no.stu, decreasing = TRUE),]
head(PestMost_stud, n = 15)

PestMost2 <- data.frame(countstudies(test, PesticideType, ActiveIngredients))
PestMost_stud <- data.frame(PestMost2)[order(PestMost2$no.stu, decreasing = TRUE),]
head(PestMost_stud, n = 15)


### SUPPL TABLE S1---------
# create a supplement table to show pesticides & their respective categories; with PollutioNameCorrected
write.csv(PestMost, "Tables/SuppTabl_PestiCategories.csv")

# How many taxa? 
length(levels(dat$TaxaGroupGSBA))
levels(dat$TaxaGroupGSBA)

# without considering unknown taxa/multi groups
length(levels(dat$TaxaGroupGSBA)) - 
  length(levels(droplevels(
    dat$TaxaGroupGSBA[dat$TaxaGroupGSBA %in% 
                        c("Macro-arthropods","Macrofauna", "Micro-arthropods", "Soil arthropods", "Soil fauna")])))

# Which taxa are the most studied?
TaxaMost <- data.frame(countstudies(dat, TaxaGroupGSBA))
TaxaMost[order(TaxaMost$no.stu, decreasing = TRUE),]

# How well are the different body size categories represented? 
plot(dat$BodySizeCat)
table(dat$BodySizeCat, dat$TaxaGroupGSBA)

countstudies(dat, BodySizeCat)


# create a supplement table to show taxa groups & their respective body size and atlas groupings
TaxaClassif <- data.frame(countstudies(dat, BodySizeCat, TaxaGroupGSBA,TaxaGroup ))


write.csv(TaxaClassif, "Tables/SuppTabl_TaxaGroupCategories.csv")


# what community metrics?
countstudies(dat, BiodivMetric)

# aggregate them into broader categories
dat$BiodivMetricBroad <- factor(ifelse(dat$BiodivMetric %in% c("Diversity indices", "Evenness indices", "Richness"), "Diversity", "Abundance & Biomass"))

countstudies(dat, BiodivMetricBroad)


# where the studies were conducted? 
summary(dat$Country)


# other metadata
countstudies(dat, ExperimentObservation)
countstudies(dat, dat$SpatialExtent.km2.) # 8 studies are lab experiments

  
  
# Pesticide inputs levels and length--------------------
summary(dat$RecommendedDose)  # raw data
summary(dat$RecommendedRate) # cleaned data


countstudies(dat, RecommendedRate)  # 31 studies at reco rate
countstudies(dat, RecommendedRate, PesticideType)
countstudies(dat, RecommendedRate, BodySizeCat)
countstudies(dat, RecommendedRate, BiodivMetric)

#### Variety of taxa groups covered by each pesticide type------------------

dat %>% group_by(PesticideType) %>% summarize(length(levels(droplevels(TaxaGroupGSBA))))

#### Short term/Long term studies------------------------------------------------

# temporal extent
countstudies(dat, TemporalExtent)
countstudies(dat, RepeatedApplication)

# repeated/single application against temporal extent
table(dat$TemporalExtent, dat$RepeatedApplication)

# we have more than 3 studies per factor level combination
countstudies(dat, TemporalExtent, RepeatedApplication)

# is there enough studies to test interaction between pesticide type and duration? yes, but some combinations have less than 3 studies for certain combination (broad spectrum, fungicides and herbicides)
countstudies(dat, TemporalExtent, PesticideType)

plot(dat$ES2~dat$TemporalExtent)


#### Multiple substances------------------------------------------------

countstudies(dat %>% filter(PesticideType=="MultiSubstance"))

# how many include insecticides?
countstudies(dat %>% filter(PesticideType=="MultiSubstance"), MultiSubstanceWithInsecticide)

# what types of pesticides?
dat %>% filter(PesticideType=="MultiSubstance") %>% group_by(PollutionNameCorrected) %>% summarize(n())

# may use with caution previous categorization
dat %>% filter(PesticideType=="MultiSubstance") %>% group_by(PollutantClass, MultiSubstanceWithInsecticide) %>% summarize(n())

countstudies(dat %>% filter(PesticideType=="MultiSubstance"), PollutantClass, MultiSubstanceWithInsecticide)
# potentially 6 studies with mixtures of herbicides; 3 studies with mix fungicides; but only 1 study with mixtures of insecticides


## Independence of covariates ---------------------------------------------------------

# what types of moderators can be included and combined into the models?

## 0.1. Pesticide type -------------------------------------------------------

# Pesticide types versus other covariates
countstudies(dat, PollutantClass1, BiodivMetricBroad)
countstudies(dat, PollutantClass1, BodySizeCat)
countstudies(dat, PollutantClass1, TaxaGroupGSBA) %>% as.data.frame()
countstudies(dat, PollutantClass1, TaxaGroupGSBA) %>% 
  group_by(PollutantClass1) %>% 
  summarize(n()) # shows that taxa range for insecticides is much broader than the other pesticide types
countstudies(dat, PollutantClass1, LongTerm2)

# other covariates, secondary
countstudies(dat, PollutantClass1, BiodivMetric)
countstudies(dat, PollutantClass1, SecondaryPollution)
countstudies(dat, PollutantClass1, SoilPH)



## 0.2. Community metric---------------------------------------------------------

# taxa breadth of biodiversity metrics: much higher for abundance: only the main taxa are represented in diversity measures. 
countstudies(dat, BiodivMetricBroad, TaxaGroupGSBA) %>% 
  group_by(BiodivMetricBroad) %>% 
  summarize(n())

# by taxa
# countstudies(dat, BiodivMetricBroad, TaxaGroupGSBA) %>% View()
countstudies(dat, BiodivMetricBroad, TaxaGroupGSBA) %>% filter(no.stu>=3)

# by pollutant: abundance have more obs of insecticides effects, while diversity has more observations of fungicide effects
countstudies(dat, BiodivMetricBroad, PollutantClass1)



## 0.3. Taxa, body size types----------------------------------------------------

# shows that cannot test interaction body size x pesticide type (only for insecticides & herbicides)
countstudies(dat, PollutantClass1, BodySizeCat) %>% filter(no.stu>=3)

# possible to test interaction biodiversity metric x body size group
countstudies(dat, BiodivMetricBroad, BodySizeCat) %>% filter(no.stu>=3)

# three wya interaction: not possible
countstudies(dat, PollutantClass1, BiodivMetricBroad, BodySizeCat) %>% filter(no.stu>=3)



## 0.4. Short term vs long term studies-------------------------------------------

# shows that possible to test interaction between type of pesticide and temporal scale
countstudies(dat, PollutantClass1, LongTerm2)

# shows that all body size categories have short, mid and long term studies
countstudies(dat, LongTerm2, BodySizeCat)

# shows that it is possible to test interaction between community metric and temporal scale
countstudies(dat, LongTerm2, BiodivMetricBroad)

# 
countstudies(dat, LongTerm2, RepeatedApplication)

#
countstudies(dat, PollutantClass1, LongTerm2, RepeatedApplication) %>% as.data.frame()

# mechanical weeding
dat[grepl("mechanical", dat$TreatmentDescription),]
dat[grepl("hand", dat$TreatmentDescription),]


## 0.5. Specific pesticides -----------------------

# can we test the effects of glyphosate and neonicotinoids? they are persistent

### 1. glyphosate -------------
summary(dat$PollutionNameCorrected)

dgly <- dat %>% filter(grepl("glyphosate", PollutionNameCorrected))

countstudies(dgly)

# 8 studies with glyphosate alone, 3 in a mixture
countstudies(dgly, PesticideType)

# only one study on microfauna
countstudies(dgly, BodySizeCat)

# mostly earthworms and micro-arthropods
countstudies(dgly, TaxaGroupGSBA)

# both repeated and single applications represented, and all time spans represented
countstudies(dgly, RepeatedApplication)
countstudies(dgly, TemporalExtent)

summary(droplevels(dgly$TreatmentDescription))

### 2. neonics -------------
summary(dat$Neonicotinoid)

dneo <- dat %>% filter(Neonicotinoid =="yes")

countstudies(dneo)

# 3 studies with neonic alone, 5 in a mixture
countstudies(dneo, PesticideType)

# only one study on microfauna
countstudies(dneo, BodySizeCat)

# mostly earthworms and micro-arthropods
countstudies(dneo, TaxaGroupGSBA)

# both repeated and single applications represented, and short and long term (but not months)
countstudies(dneo, RepeatedApplication)
countstudies(dneo, TemporalExtent)

summary(droplevels(dneo$TreatmentDescription))



## 0.6. Separate taxa groups by soft bodies/exoskeleton -----------------------

summary(dat$TaxaGroupGSBA)

dat <- dat %>% mutate(BodyTrait = factor(
  ifelse(TaxaGroupGSBA %in% c("Earthworms", "Enchytraeids", "Nematodes","Oligochaetes"), "Soft body", 
         ifelse(TaxaGroupGSBA %in% c("Macrofauna", "Soil fauna"), NA,
                "Exoskeleton"))))

countstudies(dat, BodyTrait)  # about 30 studies in each
countstudies(dat, BodyTrait, PesticideType)  # good balance across pesticide types
countstudies(dat, BodyTrait, PesticideType, BiodivMetricBroad)  # good balance across pesticide types


## References---------------------------

# list of recent studies
levels(droplevels(dat$NameOfPDF[which(grepl("2014|2015|2016|2017|2018", dat$NameOfPDF))]))



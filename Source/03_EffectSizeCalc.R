# This script computes effect sizes for each study 
# Effect sizes are Hedge's g
# Missing variances for some studies are imputed based on variances from other studies based on Koricheva2013 Handbook

## Functions ------------------------------------
library(metafor)
library(ggplot2)
library(dplyr)
library(car)


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
dat <- read.csv("Data/02PesticidesEffectsClean.csv")


## Effect sizes------------------------------------
EffectSizes <- escalc(measure = "SMD", # Hedges' d
                      m2i = Control_mean, # group 2 corresponds to the control group
                      sd2i = Control_SD,
                      n2i = Control_N,
                      
                      m1i = Treatment_mean, # group 1 is the treatment group
                      sd1i = Treatment_SD,
                      n1i = Treatment_N,
                      
                      data = dat)


# we now have the log response ratio for each case
summary(EffectSizes$yi)
# with their variance
summary(EffectSizes$vi)


## Clean the data--------------------------------
ESdat <- data.frame(dat, 
                    ES = EffectSizes$yi, 
                    VarES = EffectSizes$vi)

# Clean dataframe
class(ESdat)

# 22 NAs for variance
dataNA <- ESdat[is.na(ESdat$VarES),]
# mostly Shannons (11 obs)
summary(dataNA$Measurement)

# how many studies and observations have no variance of shannon?
countstudies(ESdat[which(is.na(ESdat$VarES) & ESdat$Measurement=="Shannon"),])



# are there shannon observation that do have variance? yes, from 9 studies, total of 25 observations
countstudies(ESdat[which(!is.na(ESdat$VarES) & ESdat$Measurement=="Shannon"),])

# how many of these only reported shannon?: n = 4
# studies with Shannon, no SD
shstu <- unique(ESdat$ID[is.na(ESdat$VarES) & ESdat$Measurement=="Shannon"])
# biodiv metrics reported
ESdat[ESdat$ID %in% shstu, c("Measurement", "ID")]
# ESdat[ESdat$ID %in% shstu,]


## Approximate variance of shannon effect sizes, based on existing data ------
VarApproxdata <- ESdat[which(!is.na(ESdat$VarES) & ESdat$Measurement=="Shannon"),]

# put control and treatment means on separate rows
VarApproxdata <- data.frame(
  ID= c(ESdat$ID, ESdat$ID),
  mean = c(ESdat$Control_mean, ESdat$Treatment_mean),
  sd = c(ESdat$Control_SD, ESdat$Treatment_SD),
  TaxaGroup = factor(as.character(ESdat$TaxaGroupGSBA, ESdat$TaxaGroupGSBA)),
  Pollut = factor(as.character(ESdat$PesticideType, ESdat$PesticideType)),
  ExpeObs = factor(as.character(ESdat$ExperimentObservation, ESdat$ExperimentObservation)),
  Measurement = factor(as.character(ESdat$Measurement, ESdat$Measurement))
)

# focus on shannons and remove SD that are NAs
VarApproxdata <- VarApproxdata[which(VarApproxdata$Measurement == "Shannon" & !is.na(VarApproxdata$sd)),]

# dataset information
nrow(VarApproxdata) # 50 observations (25 controls and 25 treatments)

plot(VarApproxdata)


# calculate CV and check homogeneity
VarApproxdata <- VarApproxdata %>% 
  dplyr::mutate(CV_forimput = sd/mean)


# descriptive stats: consider removing extreme values?
summary(VarApproxdata$CV_forimput)
hist(VarApproxdata$CV_forimput)
boxplot(VarApproxdata$CV_forimput)


# influence of case study, taxonomic group and study type?
par(mfrow=c(2,2), mar = c(2,4,1,1))
boxplot(VarApproxdata$CV_forimput~VarApproxdata$ID)
boxplot(VarApproxdata$CV_forimput~VarApproxdata$TaxaGroup, main = "Taxonomic group")
boxplot(VarApproxdata$CV_forimput~VarApproxdata$ExpeObs, main = "Study type")
boxplot(VarApproxdata$CV_forimput~VarApproxdata$Pollut, main = "Pesticide")


# relationship between mean and SD (justifies the use of CV to estimate SD that are NAs)
par(mfrow = c(1,1))
with(VarApproxdata, plot(sd ~ mean))
mod.0 <- lm(sd ~ mean, data = VarApproxdata)
summary(mod.0)

plot(VarApproxdata$sd ~VarApproxdata$mean)
abline(a = mod.0$coefficients[1], b = mod.0$coefficients[2])


## Approximate SD, but note that the model is not super predictive, we will do sensitivity analysis at the end---
ESdat$Control_SD_approx <- ifelse(ESdat$Measurement != "Shannon", ESdat$Control_SD, 
                                  ifelse(!is.na(ESdat$Control_SD), ESdat$Control_SD,
                                         # if shannon and sd = NA, sd is approx by Control_mean * average CV (sd/mean) of existing data
                                         ESdat$Control_mean*mean(VarApproxdata$CV_forimput)))

ESdat$Treatment_SD_approx <- ifelse(ESdat$Measurement != "Shannon", ESdat$Treatment_SD, 
                                    ifelse(!is.na(ESdat$Treatment_SD), ESdat$Treatment_SD,
                                           # if shannon and sd = NA, sd is approx by Treatment_mean * average CV (sd/mean) of existing data
                                           ESdat$Treatment_mean*mean(VarApproxdata$CV_forimput)))


## Effect size calculation with new SD-----------------------------
EffectSizes <- escalc(measure = "SMD", # Hedge
                      m2i = Control_mean, # group 2 corresponds to the control group
                      sd2i = Control_SD_approx,
                      n2i = Control_N,
                      
                      m1i = Treatment_mean, # group 1 is the treatment group
                      sd1i = Treatment_SD_approx,
                      n1i = Treatment_N,
                      
                      data = ESdat)



## dataframe of effect sizes and covariates
ESdat <- data.frame(ESdat, 
                    ES2 = EffectSizes$yi, 
                    VarES2 = EffectSizes$vi)

# 9 NAs for variance: for evenness indices, and when variance close to zero
dataNA2 <- ESdat[is.na(ESdat$VarES2),]
nrow(dataNA2)
(dataNA2)


# # Remove NA for varES and ES
# EffectSizes <- EffectSizes2[!is.na(EffectSizes$VarES),]
# EffectSizes <- EffectSizes[!is.na(EffectSizes$ES),]


## Save the data -------

write.csv(ESdat, "Output/EffectSizes.csv")


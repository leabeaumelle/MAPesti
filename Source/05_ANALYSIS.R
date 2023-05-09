## Meta-analysis


# This script performs the meta-analysis

## Functions--------------------------------------------
library(rio)
library(metafor)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)


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

### remove NAs
dat <- dat[!is.na(dat$ES2),]
dat <- dat[!is.na(dat$VarES2),]

# aggregate community metrics into broader categories
dat$BiodivMetricBroad <- factor(ifelse(dat$BiodivMetric %in% c("Diversity indices", "Evenness indices", "Richness"), "Diversity", "Abundance & Biomass"))


# create BodyTrait for organisms with soft body vs exoskeleton
dat <- dat %>% mutate(BodyTrait = factor(
  ifelse(TaxaGroupGSBA %in% c("Earthworms", "Enchytraeids", "Nematodes","Oligochaetes"), "Soft body", 
         ifelse(TaxaGroupGSBA %in% c("Macrofauna", "Soil fauna"), NA,
                "Exoskeleton"))))

# Data distribution
hist(dat$ES2)
dotchart(dat$ES2)


## Hypotheses 1 - General model ----------------------------
# First model investigates response of abundance-biomass and diversity-richness to different pesticides.

res1 <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType*BiodivMetricBroad, 
               random =~ 1|ID/ObsID, 
               dat = dat)

res1

# model is not significant: QM(df = 9) = 15.2782, p-val = 0.0836
anova(res1)

# test interaction with LRT
res1_interaction <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType*BiodivMetricBroad, random =~ 1|ID/ObsID, dat = dat, method = "ML")

res1_noint <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad,random =~ 1|ID/ObsID, dat = dat, method = "ML")

anova(res1_interaction, res1_noint)  # LRT 3.4300; p = 0.4886= no significant interaction

# interaction not significant, refine model
res1f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad, 
                random =~ 1|ID/ObsID,
                dat = dat)

# test main effect of pesticide type: QM(df = 4) = 8.3450, p-val = 0.0797
anova(res1f, btt=c(2:5))

# test main effect of biodiversity metric: QM(df = 1) = 3.9402, p-val = 0.0471
anova(res1f, btt=c(6))


# get estimates
res1ff <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad-1, 
                random =~ 1|ID/ObsID,
                dat = dat)
res1ff

# get marginal mean across, correcting for actual proportions of studies in different categories of moderators
# ref http://www.metafor-project.org/doku.php/tips:computing_adjusted_effects
# ref https://stat.ethz.ch/pipermail/r-sig-meta-analysis/2019-January/001396.html
colMeans(model.matrix(res1f))[-1]
predict(res1f, newmods = colMeans(model.matrix(res1f))[-1]) # is slightly lower than the overall mean from a null model. 

# get mean effect sizes across combination
coef(res1ff)[1:5]+coef(res1ff)[6]

# multisubstance ab vs div
predict(res1f, newmods = rbind(c(0,0,0,1,0),
                               c(0,0,0,1,1)))

# for fungicide ab vs diversity
predict(res1f, newmods=rbind(c(1,0,0,0,0),
                             c(1,0,0,0,1)), addx=TRUE)   

# for herbicides ab vs diversity
predict(res1f, newmods=rbind(c(0,1,0,0,0),
                             c(0,1,0,0,1)), addx=TRUE)   

# for insecticides ab vs diversity
predict(res1f, newmods=rbind(c(0,0,1,0,0),
                             c(0,0,1,0,1)), addx=TRUE)   

# for broad spectrum ab vs diversity
predict(res1f, newmods=rbind(c(0,0,0,0,0),
                             c(0,0,0,0,1)), addx=TRUE)   


# some deviation from normality detected for high values but a large number of points on the qq line
qqnorm(residuals(res1f))
qqline(residuals(res1f))


# standardized deleted residuals (recommended by Viechtbauer: https://www.wvbauer.com/lib/exe/fetch.php/articles:viechtbauer2021a.pdf)

# with rma.mv better to parallelize
# parallel::detectCores() # check no. available cores on the local machine
# library(parallel)
# sdr.res1f <- rstudent.rma.mv(res1f, cluster = dat$ID, parallel = "snow", ncpus = 2)

# qqplot of standardized deleted residuals show that most data fall close to the reference line
# qqnorm(sdr.res1f$obs$resid)
# qqline(sdr.res1f$obs$resid)

# # cook's distance
# outliers <- cooks.distance(res1f, parallel = "snow", ncpus = 2)
# plot(outliers, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance")
# # one observation has higher cook's distance: observation 204
# outliers[outliers>0.6]
# 
# dat[204,] #maybe because is a positive effect size of herbicides that removing it has strong influence on the model.
# 
# # excluding influential observation 204 yields similar results
# res1f.out <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad, 
#                     random =~ 1|ID/ObsID,
#                     dat = dat[-204,])


# ## with clusters per study ID
# outlierscluster <- cooks.distance(res1f, cluster = dat$ID, parallel = "snow", ncpus = 2)
# plot(outlierscluster, type="o", pch=19, xlab="Study", ylab="Cook's Distance", xaxt="n")
# axis(side=1, at=seq_along(outlierscluster), labels=as.numeric(names(outlierscluster)))
# # shows that three studies are influential: 
# 
# outlierscluster[outlierscluster>1]
# 
# dat[dat$ID=="1061",] # another herbicide study
# 
# # excluding influential study 1061 yields similar results
# res1f.out <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad, 
#                     random =~ 1|ID/ObsID,
#                     dat = dat[dat$ID!="1061",])



### Supplement table with model results-----
write.csv(coef(summary(res1f)), "Output/ModelResults/Model1_Overall_coefs.csv")
write.csv(data.frame(sigma  = res1f$s.names, estim =  res1f$sigma2, n.levels = res1f$s.nlevels), "Output/ModelResults/Model1_Overall_sigmas.csv")
write.csv(data.frame(QE = res1f$QE, pvalue= res1f$QEp, row.names = "Test for Residual Heterogeneity:"), "Output/ModelResults/Model1_Overall_QE.csv")



### Publication bias -------------
# funnel plot indicates could be publication bias issue
funnel(res1f, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)


# Egger test for publication bias: intercept differ from 0 at p = 0.1?
# standardized effect sizes (E*) 
dat$std.es.B <- dat$ES2/sqrt(dat$VarES2)

# standard errors (SE)
dat$se.es.B <-  sqrt(dat$VarES2)

# precision (1/SE)
dat$precision.es.B <- 1/dat$se.es.B

# Egger's regressions : E* vs. 1/SE and incoporating the main covariate
eggerreg <- rma.mv(std.es.B, VarES2, 
                   random =  ~ 1 | ID/ObsID,
                   data = dat,
                   mods =~ PesticideType*precision.es.B+BiodivMetricBroad*precision.es.B)

eggerreg  # shows that none of the intercepts differ from 0: no significant asymetry

# Make table for the supplement
suptable <- data.frame(estimate = eggerreg$b, se = eggerreg$se, pval = eggerreg$pval, CI.l = eggerreg$ci.lb, CI.u = eggerreg$ci.ub)

write.csv(suptable, "Tables/SupplTable_EggersRegression.csv")


### 1.2. Post-hoc analyses: pesticide type 4 levels---------
## Recode pesticide type according to H1.2: multiple vs single, vs single broad and single insecticides

dat$PesticideType2 <- factor(
  ifelse(dat$PesticideType =="BroadSpectrum", "BroadSpectrum", 
         ifelse(dat$PesticideType == "MultiSubstance", "MultiSubstance", 
                ifelse(dat$PesticideType == "Insecticide", "Single_Insecticide", 
                       "Single_NotInsecticide"))))
countstudies(dat, PesticideType2, BiodivMetricBroad, PesticideType)


res1.2 <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType2*BiodivMetricBroad, 
                 random =~ 1|ID/ObsID, 
                 dat = dat)

res1.2

# model is significant: QM(df = 7) = 14.0229, p-val = 0.0508
anova(res1.2)

# test interaction with LRT
res1.2_interaction <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType2*BiodivMetricBroad, random =~ 1|ID/ObsID, dat = dat, method = "ML")

res1.2_noint <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType2+BiodivMetricBroad,random =~ 1|ID/ObsID, dat = dat, method = "ML")

anova(res1.2_interaction, res1.2_noint)  # LRT 2.2993  p = 0.5126 no significant interaction

# interaction not significant, refine model
res1.2f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType2+BiodivMetricBroad, 
                  random =~ 1|ID/ObsID,
                  dat = dat)

# test main effect of pesticide type: signif. QM(df = 3) = 8.2013, p-val = 0.0420
anova(res1.2f, btt=c(2:4))

# test main effect of biodiversity metric: QM(df = 1) = 4.0111, p-val = 0.0452
anova(res1.2f, btt=c(5))


# get estimates
res1.2ff <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType2+BiodivMetricBroad-1, 
                   random =~ 1|ID/ObsID,
                   dat = dat)
res1.2ff

# get marginal mean across, correcting for actual proportions of studies in different categories of moderators
predict(res1.2f, newmods = colMeans(model.matrix(res1.2f))[-1]) # 

# get mean effect sizes across combination
coef(res1.2ff)[1:4]+coef(res1.2ff)[5]

# multisubstance ab vs div
predict(res1.2f, newmods = rbind(c(1,0,0,0),
                                 c(1,0,0,1)), addx = TRUE)

# for fungicide + herbicides ab vs diversity
predict(res1.2f, newmods=rbind(c(0,0,1,0),
                               c(0,0,1,1)), addx=TRUE)   

# for insecticides ab vs diversity
predict(res1.2f, newmods = rbind(c(0,1,0,0),
                                 c(0,1,0,1)), addx = TRUE)  

# for broad spectrum ab vs diversity
predict(res1.2f, newmods=rbind(c(0,0,0,0),
                               c(0,0,0,1)), addx=TRUE)   


### Supplement table with model results-----
write.csv(coef(summary(res1.2f)), "Output/ModelResults/Model1.2_PestiTypes_coefs.csv")
write.csv(data.frame(sigma  = res1.2f$s.names, estim =  res1.2f$sigma2, n.levels = res1.2f$s.nlevels), "Output/ModelResults/Model1.2_PestiTypes_sigmas.csv")
write.csv(data.frame(QE = res1.2f$QE, pvalue= res1.2f$QEp, row.names = "Test for Residual Heterogeneity:"), "Output/ModelResults/Model1.2_PestiTypes_QE.csv")



## Hypotheses 2 - Effect of Exposure types------------------

### 2.1.RECOMMENDED RATES-----------------------

summary(dat$RecommendedRate)
countstudies(dat, RecommendedRate)

# keep only observations at recommended rate for the analysis
datrr <- dat[which(dat$RecommendedRate == "at recommended rate"),]

countstudies(datrr)  # 31 studies
countstudies(datrr, PesticideType)  # few studies for fungicides
countstudies(datrr, BodySizeCat) # balanced body size categories
countstudies(datrr, PesticideType, BodySizeCat) # for fungicides we have only mesofauna
countstudies(datrr, PesticideType, BiodivMetricBroad) # div and ab represented although very few div obs

# other moderators
countstudies(datrr, TaxaGroupGSBA) # all main groups represented
countstudies(datrr, TemporalExtent) # ok
countstudies(datrr, RepeatedApplication) # ok


# model by pesticide type
res2.1f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType, 
               random =list(~1|ID/ObsID, ~1|BiodivMetricBroad),#~1|ID/ObsID, 
               dat = datrr)

res2.1f

# model is ns QM(df = 4) = 7.4444, p-val = 0.1142
anova(res2.1f)


# some deviation from normality detected for lowest values
qqnorm(resid(res2.1f))
qqline(resid(res2.1f))


# Calculate marginal mean effect size
colMeans(model.matrix(res2.1f))[-1] # proportions of each category, intercept left out by default
predict(res2.1f, newmods = colMeans(model.matrix(res2.1f))[-1]) #marginal mean overall, close to the main analysis



### Supplement table with model results-----
write.csv(coef(summary(res2.1f)), "Output/ModelResults/Model2.1_RecoRate_coefs.csv")
write.csv(data.frame(sigma  = res2.1f$s.names, estim =  res2.1f$sigma2, n.levels = res2.1f$s.nlevels), "Output/ModelResults/Model2.1_RecoRate_sigmas.csv")
write.csv(data.frame(QE = res2.1f$QE, pvalue= res2.1f$QEp, row.names = "Test for Residual Heterogeneity:"), "Output/ModelResults/Model2.1_RecoRate_QE.csv")


### 2.2. TEMPORAL EXTENT ------------------------------

# no long term studies for broad spectrum
countstudies(dat, PesticideType, TemporalExtent)

# subset data for analysis
dattempo <- dat %>% filter(PesticideType != "BroadSpectrum") %>% droplevels()

# count studies
countstudies(droplevels(dattempo)); countstudies(dat)

# make contingency tables of all moderators
countstudies(dattempo, PesticideType, TemporalExtent)
countstudies(dattempo, PesticideType, TemporalExtent, BiodivMetricBroad) %>% as.data.frame() # short term studies insecticides = only abundance
countstudies(dattempo, PesticideType, TemporalExtent, BodySizeCat) %>% as.data.frame() # fungicides x years= only microfauna; + microfauna absent from many factor level combinations; 
countstudies(dattempo, PesticideType,TemporalExtent, RepeatedApplication)%>% as.data.frame() # herbicide long term = only multiple applications.
countstudies(dattempo, TemporalExtent, RepeatedApplication)%>% as.data.frame() # most long term are repeated applic while short term are single applic


## model
res2.2 <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType*TemporalExtent, 
               random =list(~1|ID/ObsID, ~1|BiodivMetricBroad),#~1|ID/ObsID, 
               dat = dattempo)
res2.2

# model is ns
anova(res2.2)

# test interaction with LRT
res2.2_interaction <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType*TemporalExtent, random = list(~1|ID/ObsID, ~1|BiodivMetricBroad), dat = dattempo, method = "ML")

res2.2_noint <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+TemporalExtent,random = list(~1|ID/ObsID, ~1|BiodivMetricBroad), dat = dattempo, method = "ML")

anova(res2.2_interaction, res2.2_noint)  # LRT (df = 6) 7.8323  p = 0.2506 = no significant interaction

# interaction not significant, refine model
res2.2f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+TemporalExtent, 
                random = list(~1|ID/ObsID, ~1|BiodivMetricBroad),
                dat = dattempo)

# test main effect of tempo: QM(df = 2) = 1.2622, p-val = 0.5320, ns
anova(res2.2f, btt=c(5:6))

# test main effect of pesticide type QM(df = 3) = 3.7512, p-val = 0.2896, ns
anova(res2.2f, btt=c(2:4))

# derive mean effect sizes: multisub days and years have negative significant mean effects
predict(res2.2f, newmods = rbind(c(0,0,0,0,0),
                               c(0,0,0,1,0),
                               c(0,0,0,0,1), 
                               
                               c(1,0,0,0,0),
                               c(1,0,0,1,0),
                               c(1,0,0,0,1),
                               
                               c(0,1,0,0,0),
                               c(0,1,0,1,0),
                               c(0,1,0,0,1),
                               
                               c(0,0,1,0,0),
                               c(0,0,1,1,0),
                               c(0,0,1,0,1)),
        addx = TRUE)


# model without pesticide types: years have negative effect size overall
res2.2ff <- rma.mv(yi = ES2, V = VarES2, mods =~ TemporalExtent-1, 
                 random = list(~1|ID/ObsID, ~1|BiodivMetricBroad),
                 dat = dattempo)
res2.2ff




### Supplement table with model results-----
write.csv(coef(summary(res2.2f)), "Output/ModelResults/Model2.2_Tempo_coefs.csv")
write.csv(data.frame(sigma  = res2.2f$s.names, estim =  res2.2f$sigma2, n.levels = res2.2f$s.nlevels), "Output/ModelResults/Model2.2_Tempo_sigmas.csv")
write.csv(data.frame(QE = res2.2f$QE, pvalue= res2.2f$QEp, row.names = "Test for Residual Heterogeneity:"), "Output/ModelResults/Model2.2_Tempo_QE.csv")




## Hypotheses 3 - Effect of functional traits---------------

### 3.1. BODY SIZE ----------------

# subset of the data with body size categories
datbody <- dat %>% filter(BodySizeCat %in% c("Macrofauna", "Mesofauna", "Microfauna")) %>% droplevels()

countstudies(datbody); countstudies(dat)

# countstudies: cannot test broad spectrum-macrofauna (n = 1 study); herbicide-microfauna only has (n = 2 studies)
countstudies(datbody, BodySizeCat, PesticideType)
countstudies(datbody, BodySizeCat, BiodivMetricBroad)

# remove broad spectrum data
datbody <- datbody %>% filter(PesticideType != "BroadSpectrum")

# check sample size and independence moderators
countstudies(datbody)  # 47 studies
countstudies(datbody, PesticideType)  # ok
countstudies(datbody, BodySizeCat) # balanced body size categories
countstudies(datbody, BodySizeCat, PesticideType) # ok
countstudies(datbody, BodySizeCat, BiodivMetricBroad) # ok
countstudies(datbody, BodySizeCat, PesticideType, BiodivMetricBroad) %>% as.data.frame() # not all factor level combos are represented, some pesticides types don't have diversity observations


# other moderators
countstudies(datbody, TaxaGroupGSBA) # all main groups represented
countstudies(datbody, TemporalExtent) # ok
countstudies(datbody, RepeatedApplication) # ok

# body size cat correlated to other moderators?
countstudies(datbody, BodySizeCat, TemporalExtent) # ok
countstudies(datbody, BodySizeCat, RepeatedApplication) # ok
countstudies(datbody, BodySizeCat, RecommendedRate) # ok


# model 
res3.1 <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType*BodySizeCat, 
               random =list(~1|ID/ObsID, ~1|BiodivMetricBroad),#~1|ID/ObsID, 
               dat = datbody)

res3.1

# model is ns
anova(res3.1)

# test interaction with LRT
res3.1_interaction <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType*BodySizeCat, random = list(~1|ID/ObsID, ~1|BiodivMetricBroad), dat = datbody, method = "ML")

res3.1_noint <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BodySizeCat,random = list(~1|ID/ObsID, ~1|BiodivMetricBroad), dat = datbody, method = "ML")

anova(res3.1_interaction, res3.1_noint)  # LRT 5.4484 p = 0.4877 = no significant interaction

# interaction not significant, refine model
res3.1f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BodySizeCat, 
                random = list(~1|ID/ObsID, ~1|BiodivMetricBroad),
                dat = datbody)

# test main effect of pesticide type: QM(df = 3) = 4.6536, p-val = 0.1990
anova(res3.1f, btt=c(2:4))

# test main effect of body size: QM(df = 2) = 2.4121, p-val = 0.2994
anova(res3.1f, btt=c(5:6))

# some deviation from normality detected for high values
qqnorm(resid(res3.1f))
qqline(resid(res3.1f))


### Supplement table with model results-----
write.csv(coef(summary(res3.1f)), "Output/ModelResults/Model3.1_BodySize_coefs.csv")
write.csv(data.frame(sigma  = res3.1f$s.names, estim =  res3.1f$sigma2, n.levels = res3.1f$s.nlevels), "Output/ModelResults/Model3.1_BodySize_sigmas.csv")
write.csv(data.frame(QE = res3.1f$QE, pvalue= res3.1f$QEp, row.names = "Test for Residual Heterogeneity:"), "Output/ModelResults/Model3.1_BodySize_QE.csv")




### 3.2. EXOSKELETON-----------------------------------

# subset of the data with body traits: 5 studies removed
datsoft <- dat %>% filter(!is.na(BodyTrait)) %>% droplevels()

countstudies(datsoft); countstudies(dat)

# countstudies: enough studies to test
countstudies(datsoft, PesticideType, BodyTrait)


# check sample size and independence moderators
countstudies(datsoft)  # 47 studies
countstudies(datsoft, PesticideType)  # ok
countstudies(datsoft, BodyTrait) # balanced body size categories
countstudies(datsoft, PesticideType, BodyTrait) # ok
countstudies(datsoft, PesticideType, BodyTrait, BiodivMetricBroad) %>% as.data.frame() # not all factor level combos are represented, broad spectrum don't have have diversity observations for exoskeleton organisms

# other moderators
countstudies(datsoft, TemporalExtent) # ok
countstudies(datsoft, RepeatedApplication) # ok

# body size cat correlated to other moderators?
countstudies(datsoft, BodyTrait, TemporalExtent) # ok
countstudies(datsoft, BodyTrait, RepeatedApplication) # ok
countstudies(datsoft, BodyTrait, RecommendedRate) # ok


# model
res3.2 <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType*BodyTrait, 
               random =list(~1|ID/ObsID, ~1|BiodivMetricBroad),#~1|ID/ObsID, 
               dat = datsoft)

res3.2

# model is ns
anova(res3.2)


# test interaction with LRT
res3.2_interaction <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType*BodyTrait, random = list(~1|ID/ObsID, ~1|BiodivMetricBroad), dat = datsoft, method = "ML")

res3.2_noint <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BodyTrait,random = list(~1|ID/ObsID, ~1|BiodivMetricBroad), dat = datsoft, method = "ML")

anova(res3.2_interaction, res3.2_noint)  # LRT (df = 4) 3.2919 p = 0.5102 = no significant interaction

# interaction not significant, refine model
res3.2f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BodyTrait, 
                random = list(~1|ID/ObsID, ~1|BiodivMetricBroad),
                dat = datsoft)

# test main effect of pesticide type: significant p = 0.0447
anova(res3.2f, btt=c(2:5))

# test main effect of body trait: QM(df = 1) = 1.1116, p-val = 0.2917, ns
anova(res3.2f, btt=c(6))

# some deviation from normality detected for high values
qqnorm(resid(res3.2f))
qqline(resid(res3.2f))



### Supplement table with model results-----
write.csv(coef(summary(res3.2f)), "Output/ModelResults/Model3.2_Exoskeleton_coefs.csv")
write.csv(data.frame(sigma  = res3.2f$s.names, estim =  res3.2f$sigma2, n.levels = res3.2f$s.nlevels), "Output/ModelResults/Model3.2_Exoskeleton_sigmas.csv")
write.csv(data.frame(QE = res3.2f$QE, pvalue= res3.2f$QEp, row.names = "Test for Residual Heterogeneity:"), "Output/ModelResults/Model3.2_Exoskeleton_QE.csv")




## Supplementary analyses-----------------------------------------------------

## 1. Specific pesticides------------

### GLYPHO ------------

# filter studies applying glypho, alone or in combination
dgly <- dat %>% filter(grepl("glyphosate", PollutionNameCorrected))

countstudies(dgly)
countstudies(dgly, PesticideType)
countstudies(dgly, PesticideType, BiodivMetricBroad)
countstudies(dgly, PesticideType, TemporalExtent)
countstudies(dgly, PesticideType, RepeatedApplication)
countstudies(dgly, PesticideType, TaxaGroupGSBA)  # mostly earthworms


resgly <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType-1, 
                 random =list(~1|ID/ObsID, ~1|BiodivMetricBroad), 
                 dat = dgly)
resgly
# no difference between glypho alone or in combo with other substances (QM(df = 2) = 1.3594, p-val = 0.5068)
# effect size more negative of glypho alone, but ns. (-0.49)


###  NEONICS -------------------------------
summary(dat$Neonicotinoid) # 38 observations

countstudies(dat, PesticideType,Neonicotinoid)

datneo <- dat %>% filter(Neonicotinoid=="yes")

countstudies(datneo) # 9 studies
countstudies(datneo, PesticideType) # 5 with multisub, 6 neonics alone
countstudies(datneo, SeedCoating, PesticideType)
countstudies(datneo, TaxaGroupGSBA)
countstudies(datneo, TemporalExtent) 
countstudies(datneo, RepeatedApplication) 
countstudies(datneo, BiodivMetricBroad) 


# model
resneo <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType-1, 
                 random =list(~1|ID/ObsID, ~1|BiodivMetricBroad), 
                 dat = datneo)
resneo



## 2. Temporal effect of single insecticide applications----

# we can test temporal effect of single insecticide application only (n >=2 studies)
countstudies(dat, RepeatedApplication, PesticideType, TemporalExtent)

# subset data for analysis
datinstemp <- dat %>% filter(PesticideType == "Insecticide" & RepeatedApplication == "single") %>% droplevels()

# count studies: n = 14; 103
countstudies(droplevels(datinstemp)); countstudies(dat)

# contingency tables reveal uneven diversity representation
countstudies(datinstemp, TemporalExtent)
countstudies(datinstemp, TemporalExtent, BodySizeCat) %>% as.data.frame() # all body size represented 
countstudies(datinstemp, TemporalExtent, BiodivMetricBroad) %>% as.data.frame() # diversity observations as only present in the mid-term category

# subset data for analysis
datinstemp2 <- dat %>% filter(PesticideType == "Insecticide" & RepeatedApplication == "single" & BiodivMetricBroad == "Abundance & Biomass") %>% droplevels()

# count studies: n = 12; 89
countstudies(droplevels(datinstemp2)); countstudies(dat)

# make contingency tables of all moderators
countstudies(datinstemp2, TemporalExtent)
countstudies(datinstemp2, TemporalExtent, BodySizeCat) %>% as.data.frame() # now no more microfauna in the long term studies
countstudies(datinstemp2, TemporalExtent, TaxaGroupGSBA) %>% as.data.frame() # taxonomic breadth of mid-term studies smaller than for short- and long-term

levels(droplevels(datinstemp2$PollutionNameCorrected))  # type of insecticides 
countstudies(datinstemp2, TemporalExtent, PollutionNameCorrected)  # long term studies only address one insecticide (diflubenzuron): hard to generalize from this

resinstemp2 <- rma.mv(yi = ES2, V = VarES2, mods =~ TemporalExtent, 
                 random = ~1|ID/ObsID,
                 dat = datinstemp2)
resinstemp2

# model is not significant : QM(df = 2) = 3.0112, p-val = 0.2219
anova(resinstemp2)

# mean effect sizes: stronger negative effects on the short term
resinstemp2p <- rma.mv(yi = ES2, V = VarES2, mods =~ TemporalExtent-1, 
                      random = ~1|ID/ObsID,
                      dat = datinstemp2)
resinstemp2p


# 3. Recommended rates vs diversity/abundance-------

# data exploration following HRP comment about pop declines vs mortality at recommended rates: 
# I tested if reco rates affected div and abundance differently overall (although hard to generalize due to few studies)
# I found similar trend with stronger negative effects of reco rates on diversity than abundance, although the difference is not statistically significant.

summary(dat$RecommendedRate)


# keep only observations at recommended rate for the analysis
datrr <- dat[which(dat$RecommendedRate == "at recommended rate"),]

countstudies(datrr, BiodivMetricBroad)
countstudies(datrr, BiodivMetricBroad, PesticideType)


# overall comparison of div and ab under reco rates still shows negative response of diversity
resexplo <- rma.mv(yi = ES2, V = VarES2, mods =~ BiodivMetricBroad-1, 
                   random = ~1|ID/ObsID,
                   dat = datrr)
resexplo


# but not when filtering out broad spectrum and fungicides
resexplo <- rma.mv(yi = ES2, V = VarES2, mods =~ BiodivMetricBroad-1,
                   random = ~1|ID/ObsID,
                   dat = datrr %>%  filter(PesticideType %in% c("Herbicide", "Insecticide", "MultiSubstance")))
resexplo


# with pesticide type (but  n = 1 studies for fungi and broad spec)
resexplo <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad, 
                   random = ~1|ID/ObsID,
                   dat = datrr)
resexplo

anova(resexplo) # model is significant p : 0.0436

# main effect of pesticide type: p = 0.04, signif
anova(resexplo, btt=c(2:5))

# test main effect of diversity/ab: p = 0.13, ns
anova(resexplo, btt=c(6))

# get mean effect sizes across combination
resexplo <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad-1, 
                   random = ~1|ID/ObsID,
                   dat = datrr)

# despite ns, effect sizes for diversity are much lower
data.frame(Ab = coef(resexplo)[1:5],  # abundance
           Div = coef(resexplo)[1:5]+coef(resexplo)[6]) # diversity




# removing broad spectrum and fungicide that have few studies
resexplo2 <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad, 
                   random = ~1|ID/ObsID,
                   dat = datrr %>%  filter(PesticideType %in% c("Herbicide", "Insecticide", "MultiSubstance")))
resexplo2

anova(resexplo2) # model is merely significant p : 0.0528

# main effect of pesticide type: p = 0.04, signif
anova(resexplo2, btt=c(2:3))

# test main effect of diversity/ab: p = 0.35, ns
anova(resexplo, btt=c(4))

# get mean effect sizes across combination
resexplo2 <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad-1, 
                   random = ~1|ID/ObsID,
                   dat = datrr %>%  filter(PesticideType %in% c("Herbicide", "Insecticide", "MultiSubstance")))

# despite ns, effect sizes for diversity are much lower
data.frame(Ab = coef(resexplo2)[1:3],  # abundance
           Div = coef(resexplo2)[1:3]+coef(resexplo2)[4]) # diversity


# Splitting by pesticides when we have enough studies shows that
# insecticides and multiple substances have negative effects at recommended rates on diversity, not abundance.

# herbicide ab vs. div: mean effect sizes are ns
predict(resexplo2, newmods = rbind(c(1,0,0,0),
                                   c(1,0,0,1)))
# insecticides ab vs. div: mean effect size on div is significant
predict(resexplo2, newmods = rbind(c(0,1,0,0),
                                   c(0,1,0,1)))
# multiple substances ab vs. div: mean effect size on div is significant
predict(resexplo2, newmods = rbind(c(0,0,1,0),
                                   c(0,0,1,1)))

## 4. Compare multiple substances impacts when involve several types of pesticides versus mixtures of the same type-----------------

# data exploration

datamulti <- dat %>% 
  filter(PesticideType == "MultiSubstance") %>% droplevels()

# pesticides names
summary(datamulti$PollutionNameCorrected)

# extract pollutant names and split them in different columns
ms_pestinames <- datamulti %>% 
  select(NameOfPDF, ID, ObsID, PesticideType, PollutionNameCorrected) %>%   separate(PollutionNameCorrected, as.character(c(1:18)))

head(ms_pestinames)

# match pesticide type based on our pesticide database
# pestidatabase <- read.csv("Tables/SuppTabl_PestiCategories.csv")
# write.csv(ms_pestinames, "Output/MultiSubstance_PesticidesIdentities.csv")

# manual extraction of pesticide types for each substance in the pesticide combination
ms_database <- import("Output/MultiSubstance_PesticidesIdentitiesAndTypes.xlsx")

# keep only relevant column
ms_database <- ms_database[,c(1:3,5)]
str(ms_database)

# change factors
ms_database$ObsID <- as.factor((ms_database$ObsID))
ms_database$NameOfPDF <- as.factor(as.character(ms_database$NameOfPDF))

# bind the two datasets
datamulti2 <- left_join(datamulti, ms_database, by = c("NameOfPDF","ID","ObsID"))

# count studies and observations
countstudies(datamulti2)

countstudies(datamulti2, MultiSubType)  # n = 7; 14 cover several pesticide types, n = 14; 48 observations of multiple substances of the same pesticide type; among them, herbicides clearly dominate.

countstudies(datamulti2, BiodivMetricBroad)
countstudies(datamulti2, TemporalExtent)
countstudies(datamulti2, BodySizeCat)
countstudies(datamulti2, RepeatedApplication)
countstudies(datamulti2, RecommendedRate)

# make contingency tables to see if we can estimate meaningful mean effect sizes contrasting mono and multi types
datamulti2$MultiSubType2 <- factor(ifelse(datamulti2$MultiSubType =="MultipleTypes", "MultipleTypes", "MonoType"))

countstudies(datamulti2, MultiSubType2, BiodivMetricBroad)
countstudies(datamulti2, MultiSubType2, TemporalExtent)
countstudies(datamulti2, MultiSubType2, BodySizeCat)
countstudies(datamulti2, MultiSubType2, RepeatedApplication)
countstudies(datamulti2, MultiSubType2, RecommendedRate)  # only mono type multisubstances studies apply at rates different than reco rates, so not fully comparable

# refine recommended rates
datamulti3 <- datamulti2 %>% filter(RecommendedRate == "at recommended rate")

countstudies(datamulti3, MultiSubType2, RecommendedRate, BiodivMetricBroad)

countstudies(datamulti3, MultiSubType2, MultiSubType) # mono type involve all pesticide categories, but herbicide dominate with 5 studies.

# derive the magnitude of the trend on combined div and ab
resmultionly4 <- rma.mv(yi = ES2, V = VarES2, mods =~ MultiSubType2, 
                          random =list(~1|ID/ObsID, ~1|BiodivMetricBroad), 
                          dat = datamulti3)

resmultionly4



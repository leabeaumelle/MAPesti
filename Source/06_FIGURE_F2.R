## This script makes a figure

# Figure 2 shows the grand mean effects of different pesticide types, and the mean effect sizes at recommended rates

## Functions--------------------------------------------
library(metafor)
library(dplyr)
library(ggplot2)
library(patchwork)

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

# aggregate community metrics into broader categories
dat$BiodivMetricBroad <- factor(ifelse(dat$BiodivMetric %in% c("Diversity indices", "Evenness indices", "Richness"), "Diversity", "Abundance & Biomass"))

# create BodyTrait for organisms with soft body vs exoskeleton
dat <- dat %>% mutate(BodyTrait = factor(
  ifelse(TaxaGroupGSBA %in% c("Earthworms", "Enchytraeids", "Nematodes","Oligochaetes"), "Soft body", 
         ifelse(TaxaGroupGSBA %in% c("Macrofauna", "Soil fauna"), NA,
                "Exoskeleton"))))



## Figure 1. Revised Sept 2022--------------------------------------------------
### Panel A. Overall model and overall mean effect size-------------------------
# model pest + biodiv metric
res1f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad, 
                random =~ 1|ID/ObsID,
                dat = dat)

# sample sizes
dfsize1 <- countstudies(dat, PesticideType, BiodivMetricBroad) %>% as.data.frame()

# estimates
# matrix of model estimates (dummy coded moderators) : 
matest <- rbind(
  c(0, 0, 0, 0, 0),  # broad
  c(0, 0, 0, 0, 1),
  c(1, 0, 0, 0, 0),  # fungi
  c(1, 0, 0, 0, 1),
  c(0, 1, 0, 0, 0),  # herbi
  c(0, 1, 0, 0, 1),
  c(0, 0, 1, 0, 0),  # insec
  c(0, 0, 1, 0, 1),
  c(0, 0, 0, 1, 0),  # multi
  c(0, 0, 0, 1, 1))

# make dataframe with estimtes and CI for ggplot, and overall mean effect size
dfplot1 <- data.frame(
  Metric = c(rep(c("Abundance", "Diversity"), 5), 
             "Grand mean"),
  Pesticide = c(rep("Broad spectrum", 2), rep("Fungicide", 2), rep("Herbicide", 2), rep("Insecticide", 2), rep("Multiple substances", 2), 
                "Grand mean"),
  Estimate = c(predict(res1f, newmods = matest)$pred, 
               predict(res1f, newmods = colMeans(model.matrix(res1f))[-1])$pred),
  CI_low = c(predict(res1f, newmods = matest)$ci.lb, 
             predict(res1f, newmods = colMeans(model.matrix(res1f))[-1])$ci.lb),
  CI_up = c(predict(res1f, newmods = matest)$ci.ub, 
            predict(res1f, newmods = colMeans(model.matrix(res1f))[-1])$ci.ub),
  no.stud = c(dfsize1$no.stu,
              54),
  no.obs = c(dfsize1$no.obs, 
             294))

dfplot1$Pesticide <- factor(dfplot1$Pesticide, levels = c("Broad spectrum", "Fungicide", "Herbicide", "Insecticide", "Multiple substances", "Grand mean"))


### Panel B. Recommended rates-------------------------
# refine dataset to recommended rates
datrr <- dat[which(dat$RecommendedRate == "at recommended rate"),]

# model pest + biodiv metric: because main analysis revealed no interaction
res2f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad, 
                random =~ 1|ID/ObsID,
                dat = datrr)

# sample sizes
dfsize2 <- countstudies(datrr, PesticideType, BiodivMetricBroad) %>% as.data.frame()

# estimates
# matrix of model estimates (dummy coded moderators) : 
matest <- rbind(
  c(0, 0, 0, 0, 0),  # broad
  c(0, 0, 0, 0, 1),
  c(1, 0, 0, 0, 0),  # fungi
  c(1, 0, 0, 0, 1),
  c(0, 1, 0, 0, 0),  # herbi
  c(0, 1, 0, 0, 1),
  c(0, 0, 1, 0, 0),  # insec
  c(0, 0, 1, 0, 1),
  c(0, 0, 0, 1, 0),  # multi
  c(0, 0, 0, 1, 1))

# make dataframe with estimtes and CI for ggplot, and overall mean effect size
dfplot2 <- data.frame(
  # Metric = c(rep(c("Abundance", "Diversity"), 5)),
  Metric = c(rep(c("Abundance", "Diversity (n < 2)"), 2), rep(c("Abundance", "Diversity"), 3)),
  Pesticide = c(rep("Broad spectrum", 2), rep("Fungicide", 2), rep("Herbicide", 2), rep("Insecticide", 2), rep("Multiple substances", 2)),
  Estimate = predict(res2f, newmods = matest)$pred,
  CI_low = predict(res2f, newmods = matest)$ci.lb,
  CI_up = predict(res2f, newmods = matest)$ci.ub,
  no.stud = dfsize2$no.stu,
  no.obs = dfsize2$no.obs)


### Make plot-----------
# ggplot parameters
sizetext <- 13
sizeyaxis <- 11
sizepoint <- 3
widtherrorbar <- 0.1
sizeerrorbar <- 0.4
sizeannotation <- 2.3
sizelegend <- 8

# adds sample size to xlab
mylabs1 <- paste0("(", dfplot1$no.stud, ";" , dfplot1$no.obs, ")")

# add different color for mean overall ES
mycols1 <- c("#6B2737", "#E08E45", "black")
mypch1 <- c(rep(16, 10), 18)

F1 <- dfplot1 %>%
  ggplot(., aes(x=Pesticide, y=Estimate, group = Metric, colour = Metric))+
  geom_point(position=position_dodge(0.4), size = sizepoint, shape = mypch1)+
  geom_errorbar(position=position_dodge(0.4), 
                aes(ymin=CI_low, ymax = CI_up), width = .2)+
  geom_abline(slope = 0, intercept = 0)+
  geom_vline(xintercept = 5.5, col = "black", size = 0.2)+
  annotate("text", label = mylabs1[1], x = 1-0.33, y = dfplot1$Estimate[1], col = mycols1[1], size = sizeannotation)+
  annotate("text", label = mylabs1[2], x = 1+0.33, y = dfplot1$Estimate[2], col = mycols1[2], size = sizeannotation)+
  annotate("text", label = mylabs1[3], x = 2-0.33, y = dfplot1$Estimate[3]+0.05, col = mycols1[1], size = sizeannotation)+
  annotate("text", label = mylabs1[4], x = 2+0.33, y = dfplot1$Estimate[4], col = mycols1[2], size = sizeannotation)+
  annotate("text", label = mylabs1[5], x = 3-0.33, y = dfplot1$Estimate[5]+0.08, col = mycols1[1], size = sizeannotation)+
  annotate("text", label = mylabs1[6], x = 3+0.33, y = dfplot1$Estimate[6], col = mycols1[2], size = sizeannotation)+
  annotate("text", label = mylabs1[7], x = 4-0.33, y = dfplot1$Estimate[7], col = mycols1[1], size = sizeannotation)+
  annotate("text", label = mylabs1[8], x = 4+0.33, y = dfplot1$Estimate[8], col = mycols1[2], size = sizeannotation)+
  annotate("text", label = mylabs1[9], x = 5-0.33, y = dfplot1$Estimate[9], col = mycols1[1], size = sizeannotation)+
  annotate("text", label = mylabs1[10], x = 5+0.33, y = dfplot1$Estimate[10], col = mycols1[2], size = sizeannotation)+
  annotate("text", label = mylabs1[11], x = 6+0.33, y = dfplot1$Estimate[11], col = mycols1[3], size = sizeannotation)+
  # ggtitle("Community metric")+
  
  scale_colour_manual(values =mycols1) +
  coord_flip()+
  ylim(-1.7, 1.2)+
  theme_bw() +
  ylab("Effect size")+
  theme(legend.position = "right",
        legend.title=element_blank(),
        legend.text=element_text(size = sizelegend),
        axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = sizetext))
F1

# adds sample size to xlab
mylabs2 <- paste0("(", dfplot2$no.stud, ";" , dfplot2$no.obs, ")")

# add different color for mean overall ES
mycols2 <- c("#133C55","#0892A5", "gray")
# mypch <- c(rep(16, 10))
mypch2 <- c(rep(c(16,1),2), rep(16, 6))

F2 <- dfplot2 %>%
  ggplot(., aes(x=Pesticide, y=Estimate, group = Metric, colour = Metric))+
  geom_point(position=position_dodge(0.4), size = sizepoint, shape = mypch2)+
  geom_errorbar(position=position_dodge(0.4), 
                aes(ymin=CI_low, ymax = CI_up), width = .2)+
  geom_abline(slope = 0, intercept = 0)+
  annotate("text", label = mylabs2[1], x = 1-0.33, y = dfplot2$Estimate[1], col = mycols2[1], size = sizeannotation)+
  annotate("text", label = mylabs2[2], x = 1+0.33, y = dfplot2$Estimate[2], col = mycols2[3], size = sizeannotation)+
  annotate("text", label = mylabs2[3], x = 2-0.33, y = dfplot2$Estimate[3]+0.05, col = mycols2[1], size = sizeannotation)+
  annotate("text", label = mylabs2[4], x = 2+0.33, y = dfplot2$Estimate[4], col = mycols2[3], size = sizeannotation)+
  annotate("text", label = mylabs2[5], x = 3-0.33, y = dfplot2$Estimate[5]+0.08, col = mycols2[1], size = sizeannotation)+
  annotate("text", label = mylabs2[6], x = 3+0.33, y = dfplot2$Estimate[6], col = mycols2[2], size = sizeannotation)+
  annotate("text", label = mylabs2[7], x = 4-0.33, y = dfplot2$Estimate[7], col = mycols2[1], size = sizeannotation)+
  annotate("text", label = mylabs2[8], x = 4+0.33, y = dfplot2$Estimate[8], col = mycols2[2], size = sizeannotation)+
  annotate("text", label = mylabs2[9], x = 5-0.33, y = dfplot2$Estimate[9], col = mycols2[1], size = sizeannotation)+
  annotate("text", label = mylabs2[10], x = 5+0.33, y = dfplot2$Estimate[10], col = mycols2[2], size = sizeannotation)+
  ggtitle("Recommended rates")+
  scale_colour_manual(values =mycols2) +
  coord_flip()+
  # ylim(-1.7, 1.2)+
  ylab("Effect size")+
  theme_bw() +
  theme(
        legend.position = "right",
        legend.title=element_blank(),
        legend.text=element_text(size = sizelegend),
        axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = sizetext))
F2


F1+F2+plot_annotation(tag_levels = "A")


# save a png with high res
ggsave("Figures/Figure2.png", 
       F1+F2+plot_annotation(tag_levels = "A"), 
       width = 11.7, 
       height = 4.43)

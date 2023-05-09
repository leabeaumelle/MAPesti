## This script makes a figure

# Figure 3 shows the grand mean effects of pesticides for different study durations and different functional groups of soil invertebrates

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


# Get data for the plot ----------
## Panel A - Body size---------

# data subset: exclude studies spanning multiple body size groups, and broad spectrum studies
datbody <- dat[which(dat$BodySizeCat %in% c("Microfauna", "Mesofauna", "Macrofauna")&dat$PesticideType!="BroadSpectrum"),]

 
# mean effect sizes for the three body size categories
resbody <- rma.mv(yi = ES2, V = VarES2, mods =~ BodySizeCat-1,
                  random =list(~1|ID/ObsID, ~1|BiodivMetricBroad),
                  dat = datbody)
resbody

# sample sizes
dfsize1 <- countstudies(datbody, BodySizeCat) %>% as.data.frame()

# make dataframe with estimtes and CI for ggplot
dfplot1 <- data.frame(
  BodySize = c("Macrofauna", "Mesofauna", "Microfauna"), 
  Estimate = resbody$b, 
  CI_low = resbody$ci.lb, 
  CI_up = resbody$ci.ub, 
  no.stud = dfsize1$no.stu,
  no.obs = dfsize1$no.obs)




## Panel B - Exoskeleton ---------

# data subset: remove studies spanning multiple groups with/without exoskeletong
dataexo <- dat %>% filter(!is.na(BodyTrait))

# model exoskeleton impact
resexo <- rma.mv(yi = ES2, V = VarES2, mods =~ BodyTrait-1, 
                 random =list(~1|ID/ObsID, ~1|BiodivMetricBroad),
                 dat = dataexo)
resexo


# sample sizes
dfsize2 <- countstudies(dataexo %>% filter(!is.na(BodyTrait)), BodyTrait) %>% as.data.frame()

# make dataframe with estimtes and CI for ggplot
dfplot2 <- data.frame(
  Exoskeleton = c("presence", "absence"), 
  Estimate = resexo$b, 
  CI_low = resexo$ci.lb, 
  CI_up = resexo$ci.ub, 
  no.stud = dfsize2$no.stu,
  no.obs = dfsize2$no.obs)




## Panel C- Temporal extent ------------------------

# data subset 
dattemp <- dat %>% filter(PesticideType != "BroadSpectrum")

# model just tempo, problem is no year data for broad spectrum; dunno if want to remove it or not, if we keep it  short term have significant neg effect, if we remove it, long term effects are only sig negative. annoying
rest <- rma.mv(yi = ES2, V = VarES2, mods =~ TemporalExtent-1, 
               random =list(~1|ID/ObsID, ~1|BiodivMetricBroad),
               dat = dattemp)
rest

# sample sizes
dfsize3 <- countstudies(dattemp %>% filter(PesticideType != "BroadSpectrum"), TemporalExtent) %>% as.data.frame()

# make dataframe with estimtes and CI for ggplot
dfplot3 <- data.frame(
  TemporalExtent = c("short-term", "intermediate-term", "long-term"), 
  Estimate = rest$b, 
  CI_low = rest$ci.lb, 
  CI_up = rest$ci.ub, 
  no.stud = dfsize3$no.stu,
  no.obs = dfsize3$no.obs)





# FIGURE-----
# ggplot parameters
sizetext <- 12
sizeyaxis <- 10
sizepoint <- 3
widtherrorbar <- 0.1
sizeerrorbar <- 0.4
sizeannotation <- 2.3
sizelegend <- 8


## make panel A - Body size--------------
mylabs1 <- paste0(dfplot1$BodySize,"(", dfplot1$no.stud, ";" , dfplot1$no.obs, ")")
mycols1 <- c("#345995","#20A39E",  "#9FD356")


F3A <- dfplot1 %>%
  ggplot(., aes(x=BodySize, y=Estimate, group = BodySize, colour = BodySize))+
  geom_point(position=position_dodge(0.4), size = sizepoint)+
  geom_errorbar(position=position_dodge(0.4), 
                aes(ymin=CI_low, ymax = CI_up), width = .2)+
  geom_abline(slope = 0, intercept = 0)+
  ggtitle("Body Size")+
  scale_colour_manual(values =mycols1) +
  scale_x_discrete(breaks=waiver(), labels =mylabs1) +
  coord_flip()+
  ylim(-1, 0.5)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text=element_text(size = sizelegend),
        axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = sizetext))

F3A

## make panel B - Exoskeleton----
mylabs2 <- paste0(dfplot2$Exoskeleton,"(", dfplot2$no.stud, ";" , dfplot2$no.obs, ")")
mycols2 <- c("#2D3047","#EF6461")


F3B <- dfplot2 %>%
  ggplot(., aes(x=Exoskeleton, y=Estimate, group = Exoskeleton, colour = Exoskeleton))+
  geom_point(position=position_dodge(0.4), size = sizepoint)+
  geom_errorbar(position=position_dodge(0.4), 
                aes(ymin=CI_low, ymax = CI_up), width = .2)+
  geom_abline(slope = 0, intercept = 0)+
  ggtitle("Exoskeleton")+
  scale_colour_manual(values =mycols2) +
  scale_x_discrete(breaks=waiver(), labels =mylabs2) +
  coord_flip()+
  ylim(-1, 0.5)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text=element_text(size = sizelegend),
        axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = sizetext))



F3B

## make panel C - Temporal extent----------

mylabs3 <- paste0(dfplot3$TemporalExtent,"(", dfplot3$no.stud, ";" , dfplot3$no.obs, ")")
mycols3 <- c("black","darkgray", "lightgray")

F3C <- dfplot3 %>%
  ggplot(., aes(x=TemporalExtent, y=Estimate, group = TemporalExtent, colour = TemporalExtent))+
  geom_point(position=position_dodge(0.4), size = sizepoint)+
  geom_errorbar(position=position_dodge(0.4), 
                aes(ymin=CI_low, ymax = CI_up), width = .2)+
  geom_abline(slope = 0, intercept = 0)+
  ggtitle("Temporal extent")+
  scale_colour_manual(values =mycols3) +
  scale_x_discrete(breaks=waiver(), labels =mylabs3) +
  coord_flip()+
  ylim(-1, 0.5)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text=element_text(size = sizelegend),
        axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = sizetext))



F3C

# F3C+F3A+F3B+ plot_annotation(tag_levels = 'A')
# Width/Height ideal: <img id="plot" width="100%" height="100%" src="plot_zoom_png?width=1064&amp;height=338">

ggsave("Figures/Figure3.png", 
       F3C+F3A+F3B+ plot_annotation(tag_levels = 'A'), 
       width = 10.64, 
       height = 3.38)

# Script making supplementary figures


# data-------------
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


# Map-----------
Map <- countstudies(dat, Country)
Map[order(Map$no.obs, decreasing = TRUE),]

### Make a map ------------------------------------------------
library(maptools)
library(ggmap)
library(sf)
library(rnaturalearth)
library(countrycode)
library(viridis)
library(ggplot2)
library(reshape2)

stud <- dat %>% 
  group_by(Country) %>% 
  summarize(no.stu = length(unique(ID)))

multicountrystud <- stud[grepl(",", stud$Country),]
multicountrystud

# Manually add multiple countries studies:
addingcountries <- data.frame(rbind(c("The Netherlands", 1), c("UK", 1), c("Portugal", 1), c("Germany", 1),
                                    c("Brazil", 1), c("Portugal", 1), 
                                    c("The Netherlands", 1), c("UK", 1)))
colnames(addingcountries) <- colnames(stud)
addingcountries$no.stu <- as.numeric(addingcountries$no.stu)

stud <- rbind(as.data.frame(stud), addingcountries)

# harmonize country names
stud <- stud %>%
  mutate(iso =  countrycode(stud$Country, origin = "country.name", destination = "iso3c")) %>%
  dplyr::group_by(iso) %>%
  summarize(no.stud = sum(no.stu, na.rm = TRUE))


wm <- ne_countries(scale = 110, returnclass = "sf") %>%
  left_join(stud, by = c("adm0_a3" = "iso"))

wm <- sf::st_transform(wm,"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

mybreaks <- c(2, 4, 6, 8, 10)

map <- ggplot(wm)+
  geom_sf(aes(fill = (no.stud)))+
  scale_fill_viridis(option = "viridis", 
                     begin = 0,
                     end = 1,
                     na.value = "gray", breaks = mybreaks,
                     name = "Number of\nstudies")+
  theme_bw()+
  theme(legend.position = c(0,0),
        # legend.position = "left",
        legend.justification = c(-0.2, -0.1),
        legend.title = element_text(size=15),
        legend.text = element_text(size=14))

map

# <img id="plot" width="100%" height="100%" src="plot_zoom_png?width=997&amp;height=602">
ggsave("Figures/FigSuppl_Map.png", map, width = 9.97, height = 6.02)





# Funnel plot -----------------------------

# main model
res1f <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType+BiodivMetricBroad, 
                random =~ 1|ID/ObsID,
                dat = dat)

# funnel plot
funnel(res1f, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)

# <img id="plot" width="100%" height="100%" src="plot_zoom_png?width=880&amp;height=685">

# save a png with high res
ppi <- 300
w <- 7 # width in cm

png("Figures/FigSuppl_FunnelPlot.png",
    width=w,
    height=w*6.85/8.8,
    units = "in",
    res=ppi)

funnel(res1f, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)

dev.off()





# Glyphosate and neonicotinoids-----------
datneo <- dat %>% filter(Neonicotinoid=="yes")

resneo <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType-1, 
                 random =list(~1|ID/ObsID, ~1|BiodivMetricBroad), 
                 dat = datneo)

# sample sizes
dfsize <- countstudies(datneo, PesticideType) %>% as.data.frame()

# estimates
dfplot <- data.frame(Pesticide = c("Neonicotinoids", "Neonicotinoids + others"),
                     Estimate = c( resneo$beta), 
                     CI_low = c(resneo$ci.lb),
                     CI_up = c( resneo$ci.ub), 
                     no.stud = c(dfsize$no.stu),
                     no.obs = c(dfsize$no.obs))


# ggplot parameters
sizetext <- 13
sizeyaxis <- 9
sizepoint <- 3
widtherrorbar <- 0.1
sizeerrorbar <- 0.4

# adds sample size to xlab
myxlabs <- paste0(dfplot$Pesticide, " (", dfplot$no.stud, ";" , dfplot$no.obs, ")")

F1 <- dfplot %>%
  ggplot(., aes(x=Pesticide, y=Estimate))+
  # coord_flip()+
  geom_point( size = 3)+#, col = mycols, pch = mypch)+
  geom_errorbar(aes(ymin=CI_low, ymax = CI_up), width = .2)+#, col=mycols)+
  geom_abline(slope = 0, intercept = 0)+
  geom_vline(xintercept = 5.5, col = "black", size = 0.1)+
  scale_x_discrete(breaks=waiver(), labels =myxlabs) +
  coord_flip()+
  # ylim(-1.8,1)+
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = sizetext))
F1


# panel glyphosate
dgly <- dat %>% filter(grepl("glyphosate", PollutionNameCorrected))

resgly <- rma.mv(yi = ES2, V = VarES2, mods =~ PesticideType-1, 
                 random =list(~1|ID/ObsID, ~1|BiodivMetricBroad), 
                 dat = dgly)
resgly

# sample sizes
dfsize <- countstudies(dgly, PesticideType) %>% as.data.frame()

# estimates
dfplot <- data.frame(Pesticide = c("Glyphosate", "Glyphosate + others"),
                     Estimate = c( resgly$beta), 
                     CI_low = c(resgly$ci.lb),
                     CI_up = c( resgly$ci.ub), 
                     no.stud = c(dfsize$no.stu),
                     no.obs = c(dfsize$no.obs))


# adds sample size to xlab
myxlabs <- paste0(dfplot$Pesticide, " (", dfplot$no.stud, ";" , dfplot$no.obs, ")")

F1b <- dfplot %>%
  ggplot(., aes(x=Pesticide, y=Estimate))+
  geom_point( size = 3)+
  geom_errorbar(aes(ymin=CI_low, ymax = CI_up), width = .2)+
  geom_abline(slope = 0, intercept = 0)+
  geom_vline(xintercept = 5.5, col = "black", size = 0.1)+
  scale_x_discrete(breaks=waiver(), labels =myxlabs) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = sizetext))
F1b


# Multi substances: mono type vs multiple types ------------

datamulti <- dat %>% 
  filter(PesticideType == "MultiSubstance") %>% droplevels()

# get data of pesticide types for each substance in the pesticide combination
ms_database <- import("Output/MultiSubstance_PesticidesIdentitiesAndTypes.xlsx")
ms_database <- ms_database[,c(1:3,5)]

# change factors
ms_database$ObsID <- as.factor((ms_database$ObsID))
ms_database$NameOfPDF <- as.factor(as.character(ms_database$NameOfPDF))

# bind the two datasets
datamulti2 <- left_join(datamulti, ms_database, by = c("NameOfPDF","ID","ObsID"))

# categorize multiple substances mono vs multi
datamulti2$MultiSubType2 <- factor(ifelse(datamulti2$MultiSubType =="MultipleTypes", "MultipleTypes", "MonoType"))

# refine recommended rates
datamulti3 <- datamulti2 %>% filter(RecommendedRate == "at recommended rate")

countstudies(datamulti3, MultiSubType2, RecommendedRate, BiodivMetricBroad)


# derive the magnitude of the trend on combined div and ab
resmultionly4 <- rma.mv(yi = ES2, V = VarES2, mods =~ MultiSubType2, 
                        random =list(~1|ID/ObsID, ~1|BiodivMetricBroad), 
                        dat = datamulti3)

resmultionly4

# for MS: monotype ab vs diversity
mono <- predict(resmultionly4, newmods=0, addx=TRUE)   

# for MS: multiple types ab vs diversity
multi <- predict(resmultionly4, newmods=1, addx=TRUE)   

# sample sizes
dfsize1 <- countstudies(datamulti3, MultiSubType2) %>% as.data.frame()

# make dataframe with estimtes and CI for ggplot
dfplot1 <- data.frame(
  MultipleSubstancesOf = c("same type","several types"), 
  Estimate = c(mono$pred,multi$pred),
  CI_low = c(mono$ci.lb, multi$ci.lb),
  CI_up = c(mono$ci.ub, multi$ci.ub),
  no.stud = dfsize1$no.stu,
  no.obs = dfsize1$no.obs)

# ggplot parameters
sizetext <- 13
sizeyaxis <- 11
sizepoint <- 3
widtherrorbar <- 0.1
sizeerrorbar <- 0.4
sizeannotation <- 2.3
sizelegend <- 10

# adds sample size to xlab
myxlabs1 <- paste0(dfplot1$MultipleSubstancesOf, "(",dfplot1$no.stud, ";" , dfplot1$no.obs, ")")

F1c <- dfplot1 %>%
  ggplot(., aes(x=MultipleSubstancesOf, y=Estimate))+
  geom_point(size = sizepoint)+
  geom_errorbar(aes(ymin=CI_low, ymax = CI_up), width = .2)+
  geom_abline(slope = 0, intercept = 0)+
  geom_vline(xintercept = 5.5, col = "black", size = 0.1)+
  scale_x_discrete(breaks=waiver(), labels =myxlabs1) +
  coord_flip()+
  ggtitle("Multiple pesticides (recommended rates)")+
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = sizetext))
F1c


# Make a figure with three panels to show those results
F1c+(F1/F1b)+plot_annotation(tag_levels = "A")


# <img id="plot" width="100%" height="100%" src="plot_zoom_png?width=1056&amp;height=497">

ggsave("Figures/FigSuppl_MultiANDGlypho.png",F1c+(F1/F1b)+plot_annotation(tag_levels = "A"), width = 10.56, height = 4.97)


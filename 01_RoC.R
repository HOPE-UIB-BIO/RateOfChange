###############################################################################
### ------------ Adaptation of RATEPOL program writen in fortran ---------- ###
###############################################################################

# -----------------------------------------------------------------------------
#                           PURPOSE OF PROGRAM:
# Given a species data file of pollen data and corresponding samples ages,
# or comparable data, evaluate the rate of change in species composition
# by evaluation of dissimilarity coefficients between adjacent samples
# (large DCs being indicative of rapid change, for comparable age intervals).         
# ------------------------------------------------------------------------------

#   outline:
#  1)	Smoothing of pollen data (smooth each species using one of 5 selected smoothing method)
#  2)	(optional) Creation of time BINs (created a muster of time BINS)
#  3)	Repeat this loop (for number of randomisations):
#     a.	Sample a single age sequence from age uncertainties for all samples
#     b.	Select Units for calculation (The RoC is going to be calculated between selected Units)
#     c.	(optional) Data standardisation (subsampling pollen data in each Unit to total number of X pollen grains).
#     d.	Standardise pollen data in each Unit to proportions 
#     e.	Excluding Units and Species without any data
#     f.	Calculate DC between subsequent Units (using of the 4 selected DC method)
#     g.	Standardise the DC by age difference of Units
#  4)	Validate and summarise results from individual randomisations
#  5)	Exclude data beyond focus age
#  6)	Detect individual significant Peak points (using one of the 3 selected detection method)

# ----------------------------------------------
#                     SETUP
# ----------------------------------------------
library(tidyverse)
library(reshape2)
library(ggpubr)
library(doSNOW)
library(parallel)
library(foreach)
library(doParallel)
library(scales)
library(mgcv)
library(maps)

# ----------------------------------------------
#             LOAD DATA & FUNCTIONS
# ----------------------------------------------

tibble_Europe <- readRDS("~/DATA/input/Europe_data.13.05.20.rds")

files.sources <- list.files("~/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

# ----------------------------------------------
#               DATA EXPLORATION 
# ----------------------------------------------

# DATA based shoud be based on the criteria 
# 1) that each record need to span between ca 250-8000 years and 
# 2) samples have more than 150 grains Contain only relevant data for analysis

glimpse(tibble_Europe)

# ----------------------------------------------
#               COMPUTATION 
# ----------------------------------------------

s.time <- Sys.time()

tibble_Europe_Roc <-  tibble_Europe %>%
  mutate(., ROC = map2(filtered.counts,list_ages,
                       .f = function(.x,.y)
                         {
                         try(res <- fc_R_ratepol(
                           data.source.pollen = .x,
                           data.source.age = .y,
                           sm.type = "shep", 
                           N.points = 5, 
                           range.age.max = 500, 
                           grim.N.max = 9,
                           Working.Unit = "MW",
                           BIN.size = 500,
                           N.shifts = 5,
                           rand = 1000,
                           standardise = T, 
                           S.value = 150, 
                           DC = "chord",
                           interest.treshold = 8000,
                           Peak = "GAM",
                           Debug = F))
                         
                         if(!exists("res")){res <- "NA"}
                         return(res)
                         } ))

f.time <- Sys.time()
tot.time <- f.time - s.time
tot.time


# ----------------------------------------------
#                 SAVE RESULT 
# ----------------------------------------------

tibble_Europe_Roc <- tibble_Europe_Roc %>%
  dplyr::select(dataset.id, ROC) %>%
  filter(purrr::map(ROC, is_tibble)==T)


# saveRDS(tibble_Europe_Roc, file = "~/DATA/output/tibble_Europe_Roc200518.rds") 


# ----------------------------------------------
#               PLOT RESULTS 
# ----------------------------------------------

# loadRDS("~/DATA/output/tibble_Europe_Roc200518.rds")

tibble_Europe_work <- tibble_Europe %>%
  inner_join(tibble_Europe_Roc, by="dataset.id")


# variabe definition
age.treshold = 8000 
Roc.treshold = 2 
dataset.N = "4251"


# Perplot
tibble_Europe_work %>%
  select(dataset.id, collection.handle, long, lat, ROC) %>%
  unnest(cols = c(ROC)) %>%
  filter(AGE <= age.treshold) %>%
  ggplot(aes( y=ROC, 
              x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(age.treshold,0), ylim = c(0,Roc.treshold))+
  geom_hline(yintercept = seq(from=0,to=Roc.treshold, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age.treshold, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK==T),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BC)",
    y="Rate-of-Change score"
  )+
  facet_wrap(~ dataset.id)

ggsave("RESULTS/Europe/PerPlot.pdf",width = 50, height = 30, units= "cm", dpi= 600)


# Summary
ggarrange(
  tibble_Europe_Roc %>%
    select(dataset.id, collection.handle, long, lat, ROC) %>%
    unnest(cols = c(ROC)) %>%
    filter(AGE <= age.treshold) %>%
    ggplot(aes( x= AGE))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    coord_flip(xlim=c(age.treshold,0))+
    geom_vline(xintercept = seq(from=0,to=age.treshold, by=2000), color="gray80", size=0.1)+
    geom_rug(sides = "b")+
    geom_density(fill="gray80", color="gray50")+
    labs(
      x="Age (cal yr BC)",
      y="Density of the samples"
    ) ,
  tibble_Europe_Roc %>%
    select(dataset.id, collection.handle, long, lat, ROC) %>%
    unnest(cols = c(ROC)) %>%
    filter(AGE <= age.treshold) %>%
    ggplot(aes( y=ROC, 
                x= AGE))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    coord_flip(xlim=c(age.treshold,0), ylim = c(0,Roc.treshold))+
    geom_hline(yintercept = seq(from=0,to=Roc.treshold, by=0.5), color="gray80", size=0.1)+
    geom_vline(xintercept = seq(from=0,to=age.treshold, by=2000), color="gray80", size=0.1)+
    geom_line(aes(group=as.factor(dataset.id)), size=0.5, alpha = 1/10)+
    geom_smooth(color="purple", method = "loess",formula = y ~ x,  se=F)+
    geom_hline(yintercept = 0, color="purple", size=0.1)+
    labs(
      x="Age (cal yr BC)",
      y="Rate-of-Change score"
    )+
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y  = element_blank()
    ),
  tibble_Europe_Roc %>%
    select(dataset.id, collection.handle, long, lat, ROC) %>%
    unnest(cols = c(ROC)) %>%
    filter(AGE <= age.treshold) %>%
    filter(PEAK == T) %>%
    ggplot(aes( x= AGE))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    coord_flip(xlim=c(age.treshold,0))+
    geom_vline(xintercept = seq(from=0,to=age.treshold, by=2000), color="gray80", size=0.1)+
    geom_rug(sides = "b")+
    geom_density(fill="gray80", color="gray50")+
    labs(
      x="Age (cal yr BC)",
      y="Density of Peak-points"
    )+
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y  = element_blank()
    ),
  tibble_Europe_Roc %>%
    select(dataset.id, collection.handle, long, lat, ROC) %>%
    unnest(cols = c(ROC)) %>%
    filter(AGE <= age.treshold) %>%
    filter(PEAK == T) %>%
    ggplot(aes( y=ROC, 
                x= AGE))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    coord_flip(xlim=c(age.treshold,0), ylim = c(0,Roc.treshold))+
    geom_hline(yintercept = seq(from=0,to=Roc.treshold, by=0.5), color="gray80", size=0.1)+
    geom_vline(xintercept = seq(from=0,to=age.treshold, by=2000), color="gray80", size=0.1)+
    geom_point(color="green", size=2, shape=16, alpha=2/3)+
    geom_smooth(color="red", method = "loess", formula = y ~ x,  se=F)+
    labs(
      x="Age (cal yr BC)",
      y="Rate-of-Change Peak Points"
    )+
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y  = element_blank()
    ),
  ncol = 4, 
  widths = c(1,2,1,2),
  labels = c("A","B","C","D"))


ggsave("RESULTS/Europe/Summary.pdf",dpi= 600)



# Single plot
tibble_Europe_Roc %>%
  select(dataset.id, collection.handle, long, lat, ROC) %>%
  filter(dataset.id==dataset.N) %>%
  unnest(cols = c(ROC)) %>%
  filter(AGE <= age.treshold) %>%
  ggplot(aes( y=ROC, 
              x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim=c(age.treshold,0), ylim = c(0,Roc.treshold))+
  geom_hline(yintercept = seq(from=0,to=Roc.treshold, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age.treshold, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK==T),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BC)",
    y="Rate-of-Change score"
  )

# MAP
lat.dim <- c(min(tibble_Europe_Roc$lat),max(tibble_Europe_Roc$lat)) 
long.dim <- c(min(tibble_Europe_Roc$long),max(tibble_Europe_Roc$long))

tibble_Europe_Roc %>%
  mutate(
    N.RoC.points = select(.,ROC) %>%
      pluck(.,1) %>% 
      map_dbl(.,.f=function(x) {
        filter(x,AGE <= age.treshold) %>%
          filter(.,PEAK==T) %>%
          nrow() }),
    N.Samples = select(.,ROC) %>%
      pluck(.,1) %>% 
      map_dbl(.,.f=function(x) {
        filter(x,AGE <= age.treshold) %>%
          nrow(.) }),
    Ratio = N.RoC.points/N.Samples
  ) %>% 
  ggplot(aes(x = long, y = lat)) +
  borders(fill = "gray80", colour = "gray50") +
  coord_fixed(ylim = lat.dim, xlim = long.dim) +
  geom_point(aes(color=Ratio, size=N.RoC.points)) + 
  scale_color_gradient("Ration of Peak-points to total samples",low="black",high = "red")+
  scale_size( "Number of Peak-points")+
  labs(x = "Longitude", y = "Latitude")+
  theme_classic()

ggsave("RESULTS/Europe/RoC_map_Europe.pdf",dpi= 600)

# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())

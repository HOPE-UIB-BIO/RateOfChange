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

HOPE_data_smooth <- readRDS("~/DATA/input/HOPE_data_smooth20200730.RDS")

files.sources <- list.files("~/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

# ----------------------------------------------
#               DATA EXPLORATION 
# ----------------------------------------------

glimpse(HOPE_data_smooth)


Region_filter_tibble <- tibble(REGION = c("Europe","North America","Asia","Latin America","Africa","Oceania"),
                               max_limit =c(rep(8.5e3,2),rep(12e3,4)))

HOPE_data_work <- HOPE_data_smooth %>% 
  left_join(Region_filter_tibble, by="REGION")


# ----------------------------------------------
#               COMPUTATION 
# ----------------------------------------------

s.time <- Sys.time()

HOPE_data_work_ROC <-  HOPE_data_work %>%
  mutate(., ROC = pmap(list(filtered.counts,list_ages,max_limit),
                       .f = function(.x,.y,.z)
                         {
                         try(res <- fc_R_ratepol(
                           data.source.pollen = .x,
                           data.source.age = .y,
                           sm.type = "age.w", 
                           N.points = 5, 
                           range.age.max = 500, 
                           grim.N.max = 9,
                           Working.Unit = "MW",
                           BIN.size = 500,
                           N.shifts = 5,
                           rand = 1000,
                           standardise = T, 
                           S.value = 150, 
                           DC = "chisq",
                           interest.treshold = .z,
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

HOPE_data_work_ROC <- HOPE_data_work_ROC %>%
  dplyr::select(dataset.id, ROC) %>%
  filter(purrr::map(ROC, is_tibble)==T)


# saveRDS(HOPE_data_work_ROC, file = "~/DATA/output/HOPE_Roc20200804.RDS") 


# ----------------------------------------------
#               PLOT RESULTS 
# ----------------------------------------------

# loadRDS("~/DATA/output/tibble_Europe_Roc200518.rds")

tibble_plot_work <- HOPE_data_work %>%
  inner_join(HOPE_data_work_ROC, by="dataset.id")


# variabe definition
age.treshold = max(tibble_plot_work$max_limit) 
Roc.treshold = 2 
dataset.N = "374"


# Perplot
tibble_plot_work[1:100,] %>%
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
  geom_hline(yintercept = 0, size=0.1)+
  labs(
    x="Age (cal yr BC)",
    y="Rate-of-Change score"
  )+
  facet_wrap(~ dataset.id)

ggsave("RESULTS/Europe/PerPlot.pdf",width = 50, height = 30, units= "cm", dpi= 600)


# Single plot
tibble_plot_work %>%
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



# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())

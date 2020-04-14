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
# download.file("https://www.dropbox.com/s/3hp7rv03mkg4pjz/tibble_Europe_filtered05.03.20.RData?dl=1","~/DATA/input/tibble_Europe_filtered05.03.20.RData")

setwd("~/GITHUB/RateOfChange")

load("~/DATA/input/tibble_Europe_filtered05.03.20.RData")

files.sources <- list.files("~/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

# ----------------------------------------------
#               DATA EXPLORATION 
# ----------------------------------------------

# DATA based shoud be based on the criteria 
# 1) that each record need to span between ca 250-8000 years and 
# 2) samples have more than 150 grains Contain only relevant data for analysis

glimpse(tibble_Europe2)

# ----------------------------------------------
#               COMPUTATION 
# ----------------------------------------------

s.time <- Sys.time()

tibble_Europe_Roc <-  tibble_Europe2[-66,] %>%
  mutate(., ROC = map2(filtered.counts,list_ages,
                       .f = function(.x,.y)
                         {res <- fc_ratepol(
                           data.source.pollen = .x,
                           data.source.age = .y,
                           sm.type = "m.avg", 
                           N.points = 5, 
                           range.age.max = 500, 
                           grim.N.max = 9,
                           BIN = T,
                           BIN.size = 500,
                           Shiftbin = T,
                           N.shifts = 5,
                           rand = 1000,
                           standardise = T, 
                           S.value = 150, 
                           DC = "chisq",
                           interest.treshold = 8000,
                           Debug = F
                         )} ))

f.time <- Sys.time()
tot.time <- f.time - s.time
tot.time

# ----------------------------------------------
#                 SAVE RESULT 
# ----------------------------------------------


# saveRDS(tibble_Europe_Roc, file = "~/DATA/output/tibble_Europe_Roc200410.rds") 
# loadRDS("~/DATA/output/tibble_Europe_Roc200410.rds")

#tibble_Europe_Roc %>%
#  select(dataset.id, collection.handle, long, lat, ROC) %>%
#  unnest(cols = c(ROC)) %>%
#  write.csv(.,"~/DATA/output/results20200401.csv")


# ----------------------------------------------
#               PLOT RESULTS 
# ----------------------------------------------

fc_draw_RoC(tibble_Europe_Roc,type = "perplot", age.treshold = 8000, Roc.treshold = 3,Signif.value = "Peak.gam")
ggsave("RESULTS/Europe/PerPlot.pdf",width = 50, height = 30, units= "cm", dpi= 600)

fc_draw_RoC(tibble_Europe_Roc,type = "singleplot", dataset.N = 12, age.treshold = 8000)

fc_draw_RoC(tibble_Europe_Roc,type = "summary", age.treshold = 8000, Roc.treshold = 3, Signif.value = "Peak.gam")
ggsave("RESULTS/Europe/Summary.pdf",dpi= 600)

fc_draw_RoC(tibble_Europe_Roc,type = "map", age.treshold = 8000, Signif.value = "Peak.gam")
ggsave("RESULTS/Europe/RoC_map_Europe.pdf",dpi= 600)

# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())

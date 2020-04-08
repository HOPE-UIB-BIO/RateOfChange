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
#   1) POllen data smoothing (4 types of smoothing)
#   2) (optional) Creation of BINs
#   3) Random time sequence sampling from age-depth model
#   4) (optional) Subsetting samples in each BIN 
#   5) (optional) standardization of each sample (150 pollen)
#   4) DC calculation (4 types of DC calc.)
#   5) Age standardization of DC
#   6) (optional) Repetition of calculation for BIN shift 
#   7) Repetition of calclulation for X amount of runs
#   8) Summarisation of results from all runs 
#   9) Testing siginificance of each peak points (3 types of tests)

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

# ----------------------------------------------
#             LOAD DATA & FUNCTIONS
# ----------------------------------------------
# download.file("https://www.dropbox.com/s/3hp7rv03mkg4pjz/tibble_Europe_filtered05.03.20.RData?dl=1","~/DATA/tibble_Europe_filtered05.03.20.RData")

setwd("~/GITHUB/RateOfChange")

load("~/DATA/tibble_Europe_filtered05.03.20.RData")

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
                           sm.type = "age.w", 
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

tibble_Europe_Roc %>%
  select(dataset.id, collection.handle, long, lat, ROC) %>%
  unnest(cols = c(ROC)) %>%
  write.csv(.,"results20200401.csv")

# save.image("~/HOPE/GITHUB/RateOfChange/ENV20200403.RData")
# load("~/HOPE/GITHUB/RateOfChange/ENV20200403.RData")

# ----------------------------------------------
#               PLOT RESULTS 
# ----------------------------------------------

fc_draw_RoC(tibble_Europe_Roc,type = "perplot", age.treshold = 8000, Roc.treshold = 3,Signif.value = "Peak.gam")
ggsave("PerPlot.pdf",width = 50, height = 30, units= "cm", dpi= 600)

fc_draw_RoC(tibble_Europe_Roc,type = "singleplot", dataset.N = 12, age.treshold = 8000)

fc_draw_RoC(tibble_Europe_Roc,type = "summary", age.treshold = 8000, Roc.treshold = 3, Signif.value = "Peak.gam")
ggsave("Summary.pdf",dpi= 600)

fc_draw_RoC(tibble_Europe_Roc,type = "map", age.treshold = 8000, Signif.value = "Peak.gam")
ggsave("RoC_map_Europe.pdf",dpi= 600)

# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())

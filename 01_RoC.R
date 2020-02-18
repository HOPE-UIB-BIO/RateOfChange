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
#   1) standardization to 150 pollen
#   2) smoothing
#   3) data transformation
#   4) DC calculation
#   5) Age standardization of DC
#   6) result save

# Comments from 28.01.2020:
# standardization to 150 pollen
# We will not BIN the data and use between sample DC instead. 


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

# ----------------------------------------------
#             LOAD DATA & FUNCTIONS
# ----------------------------------------------
# download.file("https://www.dropbox.com/s/lovb5ef7o5dn9e1/tibble_Europe_filtered13.02.20.RData?dl=1","~/HOPE/Data/tibble_Europe_filtered13.02.20.RData")

setwd("~/HOPE/GITHUB/RateOfChange")
# "C:/Users/ondre/Dropbox/HOPE_data"

load("C:/Users/ondre/Dropbox/HOPE_data/tibble_Europe_filtered13.02.20.RData")

files.sources <- list.files("~/HOPE/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/HOPE/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

# ----------------------------------------------
#               DATA EXPLORATION 
# ----------------------------------------------

# DATA based shoud be based on the criteria 
# 1) that each record need to span between ca 250-8000 years and 
# 2) samples have more than 150 grains Contain only relevant data for analysis



# dataset 66, 71, 75, 126, 132 are broken 
# tibble_Europe2[c(1:65,67:70,72:74,76:125,127:131,133:N.datasets),]
# N.datasets <- nrow(tibble_Europe2)
# data.sub<-tibble_Europe2
# glimpse(data.sub)


glimpse(tibble_Europe2)

# ----------------------------------------------
#               COMPUTATION 
# ----------------------------------------------

s.time <- Sys.time()

tibble_Europe_Roc <-  tibble_Europe2 %>%
  mutate(., ROC = map2(filtered.counts,list_ages,
                       .f = function(.x,.y)
                         {res <- fc_ratepol(
                           data.source.pollen = .x,
                           data.source.age = .y,
                           rand = 9,
                           standardise = T, 
                           S.value = 150, 
                           sm.type = "grim", 
                           N.points = 5, 
                           range.age.max = 300, 
                           grim.N.max = 9,
                           DC = "chisq",
                           Debug = T
                         )} ))

f.time <- Sys.time()
tot.time <- f.time - s.time
tot.time

res.df.plot <- tibble_Europe_Roc %>%
  select(dataset.id, collection.handle, long, lat, ROC) %>%
  unnest(cols = c(ROC)) %>%
  select(.,-c(newage))

write.csv(res.df.plot,"results20202014.csv")

# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())

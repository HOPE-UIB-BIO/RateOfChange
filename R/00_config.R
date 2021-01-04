#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                     Project config
#
#                     Ondrej Mottl 
#                         2020
#
#----------------------------------------------------------#

#----------------------------------------------------------#
# 1. Load libraries and functions -----
#----------------------------------------------------------#

# delete existing workspace to start clean
rm(list = ls())

# Package version control
library(renv)
# renv::init()
# renv::snapshot(lockfile = "data/lock/revn.lock")
renv::restore(lockfile = "data/lock/revn.lock")

# libraries
library(tidyverse)
library(devtools)
library(glmmTMB)
library(parallel)
library(MuMIn)
library(emmeans)
library(performance)
library(RColorBrewer)
library(DataExplorer)


# instal RRatepol package download and attach
# devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package")

library(RRatepol)

#----------------------------------------------------------#
# 2. Load example data and custom function -----
#----------------------------------------------------------#

data_example <- RRatepol::example_data

files_sources <- list.files("R/functions/") 
sapply(paste0("R/functions/", files_sources, sep =""), source)

#----------------------------------------------------------#
# 3. Definition of variables -----
#----------------------------------------------------------#

# Number of simulated enviromental variables
N_env <-  4

# diversity of pollen taxat in simulated data
low_diversity <-  5
high_diversity <-  50

# position of the enviromental change in the sequence 
breaks_recent <-  c(2000, 3000)
breaks_late <-  c(5500, 6500)

# Number of simulated datasest of pollen data
N_rep <-  100

# template of time sequence with uneven distribution of points
time_seq <-  data_example$list_ages[[4]]$ages$age

# number of cores
n_cores <-  parallel::detectCores()

# value for beta family values
very_small_value <-  .Machine$double.eps*100

#----------------------------------------------------------#
# 4. Graphical setings  -----
#----------------------------------------------------------#

theme_set(theme_classic())

text_size <- 12

color_legen_segment <- brewer.pal(n = 3, name = 'Set2')
names(color_legen_segment) <- c("correct detection", "false positives")


color_legen_dataset_type <- brewer.pal(n = 4, name = 'Set1')
names(color_legen_dataset_type) <- 
  c("high density level_high richness",
    "high density level_low richness",
    "low_density level_high richness",
    "low_density level_low richness"
  )

color_legen_smooth <- brewer.pal(n = 5, name = 'Set3')
names(color_legen_smooth) <- c("None","M.avg","Grimm","Age.w","Shep")





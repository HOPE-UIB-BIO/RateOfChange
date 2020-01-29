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
#   5) result save

# Comments from 28.01.2020:
# standardization to 150 pollen
# We will not BIN the data and use between sample DC instead. 


# ----------------------------------------------
#                     SETUP
# ----------------------------------------------
library (tidyverse)

# ----------------------------------------------
#             LOAD DATA & FUNCTIONS
# ----------------------------------------------

load("~/HOPE/Data/European_tibble_24.01.20.Rdata")

files.sources <- list.files("~/HOPE/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/HOPE/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

# ----------------------------------------------
#               DATA MODIFICATON 
# ----------------------------------------------

data.all <- tibble_Europe_corrected

# ----------------------------------------------
#               EXPLORATION 
# ----------------------------------------------

data.sub <- data.all[1,]

glimpse(data.sub) 

data.temp <- fc_extract(data.sub)

purrr::map(data.temp, dim)

data.stand <- fc_standar(data.temp,150) 

rowSums(data.stand[,-1])



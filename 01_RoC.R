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

load("~/HOPE/Data/tibble_Europe_filtered29.01.20.Rdata")

files.sources <- list.files("~/HOPE/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/HOPE/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

# ----------------------------------------------
#               DATA EXPLORATION 
# ----------------------------------------------

# DATA based shoud be based on the criteria 
# 1) that each record need to span between ca 250-8000 years and 
# 2) samples have more than 150 grains Contain only relevant data for analysis


glimpse(tibble_Europe2)

data.sub <- tibble_Europe2[1,]

glimpse(data.sub) 

data.temp <- fc_extract(data.sub,standardise = T, S.value = 150)


# ----------------------------------------------
#               EXPLORATION 
# ----------------------------------------------

library(reshape2)

ggplot(data = reshape2::melt(data.temp$Pollen), 
       aes(y=value, 
           x=c(rep(1:nrow(data.temp$Pollen),ncol(data.temp$Pollen)) )))+
  theme_classic()+
  coord_flip(ylim = c(0,1))+
  geom_point(alpha=1/5)+
  geom_line(group=1)+
  facet_wrap(~variable)


test<- fc_smooth(data.source = data.temp, sm.type = "m.avg", N.points = 5)
test<- fc_smooth(data.source = data.temp, sm.type = "grim", N.points = 5,grim.N.max = 7, range.age.max = 300)
test<- fc_smooth(data.source = data.temp, sm.type = "age.w", N.points = 5, range.age.max = 300)
test<- fc_smooth(data.source = data.temp, sm.type = "shep")

ggplot(data = reshape2::melt(test$Pollen), 
       aes(y=value, 
           x=c(rep(1:nrow(test$Pollen),ncol(test$Pollen)) )))+
  theme_classic()+
  coord_flip(ylim = c(0,1))+
  geom_point(alpha=1/5)+
  geom_line(group=1)+
  facet_wrap(~variable)

# ----------------------------------------------
#               RESULT SAVE 
# ----------------------------------------------


# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())

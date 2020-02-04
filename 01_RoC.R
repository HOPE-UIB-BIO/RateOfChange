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
library (tidyverse)
library(reshape2)
library(mvpart)

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

# ----------------------------------------------
#               COMPUTATION 
# ----------------------------------------------

N.datasets <- nrow(tibble_Europe2)

data.sub<-tibble_Europe2[c(1:N.datasets),]

list.res <- vector("list",length=nrow(data.sub))

s.time <- Sys.time()
for (i in 1:nrow(data.sub)) 
{
  list.res[[i]] <- fc_ratepol(data.sub[i,],
                              standardise = T, 
                              S.value = 150, 
                              sm.type = "grim", 
                              N.points = 5, 
                              range.age.max = 300, 
                              grim.N.max = 9,
                              DC = "chisq",
                              result = "small")
}
f.time <- Sys.time()
tot.time <- f.time - s.time
tot.time

# Why purrr does not work???

# ----------------------------------------------
#             RESULT VISUALISATION 
# ----------------------------------------------

# plot creation 
res.df.plot <- data.frame(matrix(ncol = 3, nrow = 0))
names(res.df.plot) <- c("Age","Roc","ID")

for (k in 1:length(list.res))
{
  list.res[[k]]$data$ID <- rep(list.res[[k]]$ID,nrow(list.res[[k]]$data))
  res.df.plot <- rbind(res.df.plot,list.res[[k]]$data)
}

res.df.plot[res.df.plot$RoC==max(res.df.plot$RoC),]
summary(as.factor(res.df.plot$ID))



RoC_summary <- res.df.plot[is.infinite(res.df.plot$RoC)==F,] %>%
  ggplot(aes( y=RoC, 
              x= Age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip()+
  #geom_point(alpha=1/5)+
  geom_line(aes(group=as.factor(ID)),alpha=1/10)+
  geom_smooth(color="blue", method = "loess", se=F)
RoC_summary

ggsave("RoC_summary.pdf")

# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())

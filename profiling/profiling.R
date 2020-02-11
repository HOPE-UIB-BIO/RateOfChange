# setup
library(tidyverse)
library(reshape2)
library(ggpubr)
library(doSNOW)
library(parallel)
library(foreach)
library(doParallel)

library(profvis)
library(microbenchmark)

load("~/HOPE/Data/tibble_Europe_filtered29.01.20.Rdata")
files.sources <- list.files("~/HOPE/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/HOPE/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

N.datasets <- nrow(tibble_Europe2)
data.sub<-tibble_Europe2[c(1:65,67:70,72:74,76:125,127:131,133:N.datasets),]


# test for difefrent setting 

DF.performance <- data.frame(matrix(nrow = 16, ncol=5))
names(DF.performance) <- c("smooth","DC","user","system","elapsed")
DF.performance$smooth <- c(rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4))
DF.performance$DC <- c(rep(c("euc","euc.sd","chord","chisq"),4))

dataset.N <- 3

for(i in 1:nrow(DF.performance))
{
  a<- system.time(fc_ratepol(data.sub[dataset.N,],
                             rand = 9,
                             standardise = T, 
                             S.value = 150, 
                             sm.type = DF.performance$smooth[i], 
                             N.points = 5, 
                             range.age.max = 300, 
                             grim.N.max = 9,
                             DC = DF.performance$DC[i],
                             Debug = F))
  
  DF.performance$user[i] <- a[1]
  DF.performance$system[i] <- a[2]
  DF.performance$elapsed[i] <- a[3]
}

DF.performance %>%
  ggplot(aes(y=elapsed, x=smooth))+
  geom_bar(aes(fill=DC),stat="identity", position = "dodge", color="black")+
  ggtitle(paste("ID",data.sub$site.id[[dataset.N]],",N samples",nrow(data.sub$filtered.counts[[dataset.N]])))+
  theme_classic()

ggsave("ComputationTime2.pdf")

DF.performance[order(DF.performance$elapsed),]


# test of individual code parts

data.source <- tibble_Europe2[2,]
rand = 9
standardise = T
S.value = 150
sm.type = "shep" 
N.points = 5
range.age.max = 300
grim.N.max = 9
DC = "euc.sd"
Debug = F
dataset.ID <- data.source$dataset.id


profvis({
  
  data.work <- fc_extract(data.source, Debug=Debug) 
  
  data.sd <- data.work
  
  # ----------------------------------------------
  #             DATA STANDARFISATION
  # ----------------------------------------------
  # standardisation of pollen data to X(S.value) number of pollen grains 
  if(standardise==T) # 
  {
    data.sd$Pollen <- fc_standar(data.work$Pollen, S.value, Debug=Debug)
    
    if(any(rowSums(data.sd$Pollen)!=S.value))
      stop("standardisation was unsuccesfull")
  }
  
  # data check with proportioning
  data.sd <- fc_check(data.sd, proportion = T, Debug=Debug)
  
  # ----------------------------------------------
  #               DATA SMOOTHING
  # ----------------------------------------------
  # smooth pollen data by selected smoothing type
  data.smooth <- fc_smooth(data.sd, 
                           sm.type = sm.type, 
                           N.points = N.points,
                           grim.N.max = grim.N.max, 
                           range.age.max = range.age.max,
                           Debug=Debug)
  
  #data check (with proportioning ???)
  data.smooth <- fc_check(data.smooth, proportion = T, Debug=Debug)
  
  # ----------------------------------------------
  #               DC CALCULATION
  # ----------------------------------------------
  # calculate DC for each sample
  DC.res <- fc_calDC(data.smooth,DC=DC, Debug=Debug)
  
  # ----------------------------------------------
  #             AGE STANDARDISATION
  # ----------------------------------------------
  
  sample.size.work <- data.smooth$Dim.val[2]-1 
  
  age.diff <- vector(mode = "numeric", length = sample.size.work )
  for (i in 1:sample.size.work)
  {
    age.diff[i] <- abs(data.smooth$Age$newage[i+1]-data.smooth$Age$newage[i]) 
    # temporary fix for errors in age data where age difference between samples is 0
    if(age.diff[i]==0)
    {age.diff[i]<-1}
  }
  
  
  if (Debug ==T)
  {
    print("-")
    print(paste("The time standardisation unit (TSU) is",round(mean(age.diff),2)))  
  }
  
  
  DC.res.s <- vector(mode = "numeric", length = sample.size.work)
  for (j in 1:sample.size.work)
  {
    DC.res.s[j] <- DC.res[j]*mean(age.diff)/age.diff[j]
  }
  
  # ----------------------------------------------
  #         RESULT OF SINGLE RAND RUN
  # ----------------------------------------------
  
  data.result <- data.frame(Age=data.smooth$Age$age[1:sample.size.work], RoC=DC.res.s)
  row.names(data.result) <- row.names(data.smooth$Pollen)[1:sample.size.work]
  
  data.result.temp <- as.data.frame(list(ID=1,DF=data.result)) 
  
}, interval = 0.001)




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
# download.file("https://www.dropbox.com/s/lovb5ef7o5dn9e1/tibble_Europe_filtered13.02.20.RData?dl=1","~/HOPE/Data/tibble_Europe_filtered13.02.20.RData")

setwd("~/HOPE/GITHUB/RateOfChange")
# "C:/Users/ondre/Dropbox/HOPE_data"

load("C:/Users/ondre/Dropbox/HOPE_data/tibble_Europe_filtered18.02.20.RData")

files.sources <- list.files("~/HOPE/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/HOPE/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

# test for difefrent setting 
DF.performance <- data.frame(matrix(nrow = 16, ncol=5))
names(DF.performance) <- c("smooth","DC","user","system","elapsed")
performance.smooth <- c(rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4))
performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),4))
DF.performance$smooth <- performance.smooth
DF.performance$DC <- performance.DC

dataset.N <- 2

# 3 & 187

for(i in 1:nrow(DF.performance))
{
  a<- system.time(fc_ratepol( data.source.pollen =  tibble_Europe2$filtered.counts[[dataset.N]],
                              data.source.age = tibble_Europe2$list_ages[[dataset.N]],
                              sm.type = DF.performance$smooth[i] ,
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              BIN = T,
                              BIN.size = 500,
                              Shiftbin = T,
                              N.shifts = 10,
                              rand = 10,
                              standardise = T, 
                              S.value = 150 ,
                              DC = DF.performance$DC[i],
                              interest.treshold = 8000,
                              Debug = F))
  
  DF.performance$user[i] <- a[1]
  DF.performance$system[i] <- a[2]
  DF.performance$elapsed[i] <- a[3]
}

DF.performance %>%
  ggplot(aes(y=elapsed, x=smooth))+
  geom_bar(aes(fill=DC),stat="identity", position = "dodge", color="gray50")+
  ggtitle(paste("ID",tibble_Europe2$dataset.id[[dataset.N]],
                ",N samples",nrow(tibble_Europe2$filtered.counts[[dataset.N]]),
                "N.randomisation",rand))+
  theme_classic()

ggsave(paste0("Performancetest_RUN_",dataset.N,".pdf"))

DF.performance[order(DF.performance$elapsed),]


# comparison of result between diferent settings

performance.list <- vector("list",length = 16)

for(i in 1:length(performance.list))
{
  a<- fc_ratepol( data.source.pollen =  tibble_Europe2$filtered.counts[[dataset.N]],
                  data.source.age = tibble_Europe2$list_ages[[dataset.N]],
                  sm.type = DF.performance$smooth[i] ,
                  N.points = 5,
                  range.age.max = 500, 
                  grim.N.max = 9,
                  BIN = T,
                  BIN.size = 500,
                  Shiftbin = T,
                  N.shifts = 5,
                  rand = 100,
                  standardise = T, 
                  S.value = 150 ,
                  DC = DF.performance$DC[i],
                  interest.treshold = 8000,
                  Debug = F)
  
  performance.list[[i]] <- as.data.frame(a)
}

performance.list.plot <- vector("list",length = 16)

for(i in 1:length(performance.list))
{
  data.temp <- performance.list[[i]]
  
  roc.max <- max(data.temp$RUN.RoC)*1.5
  age.max <- max(data.temp$RUN.Age.Pos)
  
  performance.list.plot[[i]]<-ggplot(data=data.temp, aes(y=RUN.RoC, 
                                             x= RUN.Age.Pos))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    coord_flip(xlim=c(0,age.max), ylim = c(0,roc.max))+
    geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), alpha=1/5, color="gray")+
    geom_line(alpha=1, size=1)+
    geom_line(data=data.frame(RUN.RoC = predict.gam(gam(RUN.RoC~s(RUN.Age.Pos), data = data.temp)),
                              RUN.Age.Pos = data.temp$RUN.Age.Pos),
              color="blue", size=1)+
    geom_point(color="black", size=1)+
    geom_point(data = data.temp[data.temp$soft.Peak==T,],color="yellow", size=1)+
    geom_point(data = data.temp[data.temp$Peak==T,],color="orange", size=2)+
    geom_point(data = data.temp[data.temp$Peak.gam==T,],color="red", size=3)+
    geom_hline(yintercept = 0, color="red")+
    xlab("Age")+ylab("Rate of Change")+
    ggtitle(paste(performance.smooth[i],"+",performance.DC[i]))
  #assign(paste0("p",i),p.temp)
}

p.all <- ggarrange(performance.list.plot[[1]],performance.list.plot[[2]],performance.list.plot[[3]],performance.list.plot[[4]],
                   performance.list.plot[[5]],performance.list.plot[[6]],performance.list.plot[[7]],performance.list.plot[[8]],
                   performance.list.plot[[9]],performance.list.plot[[10]],performance.list.plot[[11]],performance.list.plot[[12]],
                   performance.list.plot[[13]],performance.list.plot[[14]],performance.list.plot[[15]],performance.list.plot[[16]],
                   nrow = 4, ncol = 4)

p.all

ggsave("Comparison.pdf", p.all,  width = 50, height = 30, units= "cm", dpi= 600)

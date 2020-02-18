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

# test for difefrent setting 
DF.performance <- data.frame(matrix(nrow = 16, ncol=5))
names(DF.performance) <- c("smooth","DC","user","system","elapsed")
performance.smooth <- c(rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4))
performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),4))
DF.performance$smooth <- performance.smooth
DF.performance$DC <- performance.DC

dataset.N <- 128
rand = 9

# 3 & 187

for(i in 1:nrow(DF.performance))
{
  a<- system.time(fc_ratepol(data.source.pollen = tibble_Europe2$filtered.counts[[dataset.N]],
                             data.source.age =  tibble_Europe2$list_ages[[dataset.N]],
                             rand = rand,
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
  ggtitle(paste("ID",tibble_Europe_Roc$dataset.id[[dataset.N]],
                ",N samples",nrow(tibble_Europe_Roc$filtered.counts[[dataset.N]]),
                "N.randomisation",rand))+
  theme_classic()

ggsave("ComputationTime.pdf")

DF.performance[order(DF.performance$elapsed),]


# comparison of result between diferent settings

performance.list <- vector("list",length = 16)

for(i in 1:length(performance.list))
{
  a<- fc_ratepol(data.source.pollen = tibble_Europe2$filtered.counts[[dataset.N]],
                data.source.age =  tibble_Europe2$list_ages[[dataset.N]],
                             rand = 9,
                             standardise = T, 
                             S.value = 150, 
                             sm.type = performance.smooth[i], 
                             N.points = 5, 
                             range.age.max = 300, 
                             grim.N.max = 9,
                             DC = performance.DC[i],
                             Debug = F)
  
  performance.list[[i]] <- as.data.frame(a)
}

performance.list.plot <- vector("list",length = 16)

for(i in 1:length(performance.list))
{
  data.temp <- performance.list[[i]]
  
  roc.max <- max(data.temp$RoC.95q)
  age.max <- max(data.temp$age)
  
  performance.list.plot[[i]]<-ggplot(data=data.temp, aes(y=RoC.median, 
                                             x= age))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    coord_flip(xlim=c(0,age.max), ylim = c(0,roc.max))+
    geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), alpha=1/5)+
    geom_line(alpha=1, size=2)+
    geom_point(data = data.temp[data.temp$Data.Peak==T,],color="blue", alpha=1, size=3)+
    geom_hline(yintercept = median(data.temp$Data.RoC.median), color="blue")+
    geom_hline(yintercept = 0, color="red")+
    xlab("Age")+ylab("Rate of Change")+
    ggtitle(paste(performance.smooth[i],"+",performance.DC[i]))
  #assign(paste0("p",i),p.temp)
}

performance.list.plot[[2]]<- performance.list.plot[[2]]+coord_flip(xlim=c(0,15000), ylim = c(0,50))
performance.list.plot[[6]]<- performance.list.plot[[6]]+coord_flip(xlim=c(0,15000), ylim = c(0,50)) 
performance.list.plot[[10]]<- performance.list.plot[[10]]+coord_flip(xlim=c(0,15000), ylim = c(0,50))
performance.list.plot[[14]]<- performance.list.plot[[14]]+coord_flip(xlim=c(0,15000), ylim = c(0,50))

p.all <- ggarrange(performance.list.plot[[1]],performance.list.plot[[2]],performance.list.plot[[3]],performance.list.plot[[4]],
                   performance.list.plot[[5]],performance.list.plot[[6]],performance.list.plot[[7]],performance.list.plot[[8]],
                   performance.list.plot[[9]],performance.list.plot[[10]],performance.list.plot[[11]],performance.list.plot[[12]],
                   performance.list.plot[[13]],performance.list.plot[[14]],performance.list.plot[[15]],performance.list.plot[[16]],
                   nrow = 4, ncol = 4)

p.all

ggsave("Comparison.pdf", p.all,  width = 50, height = 30, units= "cm", dpi= 600)

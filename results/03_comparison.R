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


performance.smooth <- c(rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4))
performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),4))
performance.list <- vector("list",length = 16)

dataset.N <- 1

for(i in 1:length(performance.list))
{
  a<- fc_ratepol(data.sub[dataset.N,],
                             rand = 99,
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

for(i in 1:length(performance.list))
{
  data.temp <- performance.list[[i]]
  p.temp<-ggplot(data=data.temp, aes(y=Data.RoC.mean, 
                                             x= Data.DF.Age))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    coord_flip(xlim=c(0,15000), ylim = c(0,0.7))+
    geom_ribbon(aes(ymin=Data.RoC.05q, ymax=Data.RoC.95q), alpha=1/2)+
    geom_line(aes(group=as.factor(ID)),alpha=1, size=1)+
    geom_point(data = data.temp[data.temp$Data.Peak==T,],color="blue", alpha=1, size=3)+
    geom_hline(yintercept = median(data.temp$Data.RoC.mean), color="blue")+
    geom_hline(yintercept = 0, color="red")+
    xlab("Age")+ylab("Rate of Change")+
    ggtitle(paste(performance.smooth[i],"+",performance.DC[i]))
  
  assign(paste0("p",i),p.temp)
}

p2<- p2+coord_flip(xlim=c(0,15000), ylim = c(0,50))
p6<- p6+coord_flip(xlim=c(0,15000), ylim = c(0,50)) 
p10<- p10+coord_flip(xlim=c(0,15000), ylim = c(0,50))
p14<- p14+coord_flip(xlim=c(0,15000), ylim = c(0,50))

p.all <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16, nrow = 4, ncol = 4)

p.all

ggsave("Comparison.pdf", p.all,  width = 50, height = 30, units= "cm", dpi= 600)

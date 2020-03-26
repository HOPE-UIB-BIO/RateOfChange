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
# download.file("https://www.dropbox.com/s/lovb5ef7o5dn9e1/tibble_Europe_filtered13.02.20.RData?dl=1","~/HOPE/Data/tibble_Europe_filtered13.02.20.RData")

setwd("~/HOPE/GITHUB/RateOfChange")
# "C:/Users/ondre/Dropbox/HOPE_data"

load("C:/Users/ondre/Dropbox/HOPE_data/tibble_Europe_filtered05.03.20.RData")

files.sources <- list.files("~/HOPE/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/HOPE/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)


glimpse(tibble_Europe2)

plot.pollen <- function (data, sm.type, N.taxa, interest.treshold)
{
  Common.list <- data$filtered.counts %>%
    colSums() %>% 
    sort(decreasing = T) %>%
    .subset(.,1:N.taxa) %>%
    names() %>%
    sub("/",".",.)
  
  data.ext <-  fc_extract(data$filtered.counts,
                          data$list_ages) %>%
    fc_smooth(.,sm.type = sm.type,
              N.points = 5,
              grim.N.max = 9,
              range.age.max = 500) %>%
    fc_check(.,proportion = T)
  
  plot.p <- data.ext$Pollen %>%
    select(Common.list) %>%
    rownames_to_column() %>%
    pivot_longer(., cols = c(Common.list)) %>%
    rename(sample.id = rowname) %>%
    inner_join(.,data.ext$Age, by="sample.id")  %>%
    ggplot(aes( y=value, 
                x= age))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
                color="gray20", alpha=1/5, size=0.1)+
    xlab("Age")+ylab("Pollen (%)")+
    coord_flip(xlim=c(0,interest.treshold), ylim = c(0,1))+
    ggtitle(paste("smoothing",sm.type))  
  return (plot.p)
  
}


plot.comparison <- function(data, BIN, BIN.size, Shiftbin, N.shifts, rand, interest.treshold)
{
  performance.list.plot <- vector("list",length = 20);
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  
  
  for(i in 1:length(performance.list.plot))
  {
    data.temp<- fc_ratepol( data.source.pollen =  data$filtered.counts,
                            data.source.age = data$list_ages,
                            sm.type = performance.smooth[i],
                            N.points = 5,
                            range.age.max = 500, 
                            grim.N.max = 9,
                            BIN = BIN,
                            BIN.size = BIN.size,
                            Shiftbin = Shiftbin,
                            N.shifts = N.shifts,
                            rand = rand,
                            standardise = F, 
                            S.value = 150 ,
                            DC = performance.DC[i],
                            interest.treshold = interest.treshold,
                            Debug = F) %>%
      as.data.frame();
    
    
    roc.max <- max(data.temp$RUN.RoC)*1.5;
    age.max <- max(data.temp$RUN.Age.Pos);
    
    performance.list.plot[[i]]<-ggplot(data=data.temp, 
                                       aes(y=RUN.RoC, 
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
      geom_point(data = data.temp[data.temp$Peak.treshold==T,],color="yellow", size=1)+
      geom_point(data = data.temp[data.temp$Peak.treshold.95==T,],color="orange", size=2)+
      geom_point(data = data.temp[data.temp$Peak.gam==T,],color="red", size=3)+
      geom_point(data = data.temp[data.temp$Peak.SNI==T,],color="purple", size=3)+
      geom_hline(yintercept = 0, color="red")+
      geom_hline(yintercept = median(data.temp$RUN.RoC), color="green")+
      xlab("Age")+ylab("Rate of Change")+
      ggtitle(paste(performance.smooth[i],"+",performance.DC[i]))
  }
  
  p.all <- ggarrange(performance.list.plot[[1]],performance.list.plot[[2]],performance.list.plot[[3]],performance.list.plot[[4]],
                     performance.list.plot[[5]],performance.list.plot[[6]],performance.list.plot[[7]],performance.list.plot[[8]],
                     performance.list.plot[[9]],performance.list.plot[[10]],performance.list.plot[[11]],performance.list.plot[[12]],
                     performance.list.plot[[13]],performance.list.plot[[14]],performance.list.plot[[15]],performance.list.plot[[16]],
                     performance.list.plot[[17]],performance.list.plot[[18]],performance.list.plot[[19]],performance.list.plot[[20]],
                     nrow = 5, ncol = 4)
  
  return(p.all)
}


# -----------------------------------------
#                 SIMULATION
# -----------------------------------------

#tibble_Europe2$list_ages[[2]]$ages$age
data.sim <-  fc_random_data(time = seq(from=0, to=10e3, by=100),
                            nforc = 4, 
                            nprox = 10,
                            manual.edit = T,
                            breaks=c(0,2000, 4000,6000),
                            jitter = T)

plot.pollen.sim <- ggarrange(
  plot.pollen(data.sim,"none",10,8000),
  plot.pollen(data.sim,"m.avg",10,8000),
  plot.pollen(data.sim,"age.w",10,8000),
  plot.pollen(data.sim,"grim",10,8000),
  plot.pollen(data.sim,"shep",10,8000),
  ncol=5, nrow = 1, common.legend = T, legend = "right"
)

plot.pollen.sim

dataset.sim.comparison.sample <- plot.comparison(data.sim,
                                                   BIN = F, 
                                                   Shiftbin = F, 
                                                   rand = 1, 
                                                   interest.treshold =  8000)
dataset.sim.comparison.BIN <- plot.comparison(data.sim,
                                              BIN = T,
                                              BIN.size = 500,
                                              Shiftbin = F, 
                                              rand = 1, 
                                              interest.treshold =  8000)
dataset.sim.comparison.BIN.shift <- plot.comparison(data.sim,
                                                 BIN = F, 
                                                 Shiftbin = F, 
                                                 rand = 1, 
                                                 interest.treshold =  8000)

signif.sim.comparison.sample <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                    nforc=4,
                                                    mean=100, 
                                                    sdev=.15, 
                                                    nprox=10, 
                                                    var=20,
                                                    range=15,
                                                    manual.edit = T,
                                                    breaks=c(2000,3000),
                                                    jitter = T,
                                                    BIN=F,
                                                    BIN.size=500, 
                                                    Shiftbin=F,
                                                    N.shifts=5,
                                                    rand.sets=10,
                                                    interest.treshold=8000)

signif.sim.comparison.BIN.shift <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                    nforc=4,
                                                    mean=100, 
                                                    sdev=.15, 
                                                    nprox=10, 
                                                    var=20,
                                                    range=15,
                                                    manual.edit = T,
                                                    breaks=c(2000,3000),
                                                    jitter = T,
                                                    BIN=T,
                                                    BIN.size=500, 
                                                    Shiftbin=T,
                                                    N.shifts=5,
                                                    rand.sets=10,
                                                    interest.treshold=8000)


signif.sim.comparison.BIN <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                 nforc=4,
                                                 mean=100, 
                                                 sdev=.15, 
                                                 nprox=10, 
                                                 var=20,
                                                 range=15,
                                                 manual.edit = T,
                                                 breaks=c(2000,3000),
                                                 jitter = T,
                                                 BIN=T,
                                                 BIN.size=500, 
                                                 Shiftbin=F,
                                                 N.shifts=5,
                                                 rand.sets=10,
                                                 interest.treshold=8000)


# -----------------------------------------
#         SIMULATION with template
# -----------------------------------------

# uneven distribution of samples
tibble_Europe2$list_ages[[2]]$ages %>%
  filter(age <8000) %>%
  nrow

tibble_Europe2$list_ages[[2]]$ages %>%
  filter(age <8000) %>%
  ggplot(aes(x=age))+
  geom_density(fill="gray80",color="gray30")+
  coord_cartesian(xlim=c(0,8000))

data.sim.template <-  fc_random_data(time = tibble_Europe2$list_ages[[2]]$ages$age,
                            nforc = 4, 
                            nprox = 10,
                            manual.edit = T,
                            breaks=c(0,2000, 4000,6000),
                            jitter = T)

plot.pollen.sim.template <- ggarrange(
  plot.pollen(data.sim,"none",10,8000),
  plot.pollen(data.sim,"m.avg",10,8000),
  plot.pollen(data.sim,"age.w",10,8000),
  plot.pollen(data.sim,"grim",10,8000),
  plot.pollen(data.sim,"shep",10,8000),
  ncol=5, nrow = 1, common.legend = T, legend = "right"
)

plot.pollen.sim.template

dataset.sim.template.comparison.sample <- plot.comparison(data.sim.template,
                                                 BIN = F, 
                                                 Shiftbin = F, 
                                                 rand = 1, 
                                                 interest.treshold =  8000)
dataset.sim.template.comparison.BIN <- plot.comparison(data.sim.template,
                                              BIN = T,
                                              BIN.size = 500,
                                              Shiftbin = F, 
                                              rand = 1, 
                                              interest.treshold =  8000)
dataset.sim.template.comparison.BIN.shift <- plot.comparison(data.sim.template,
                                                    BIN = F, 
                                                    Shiftbin = F, 
                                                    rand = 1, 
                                                    interest.treshold =  8000)



# -----------------------------------------
#                 25318
# -----------------------------------------



dataset.25318 <- tibble_Europe2 %>%
  filter(dataset.id=="25318")


# POllen graph

plot.pollen.25318 <- ggarrange(
  plot.pollen(dataset.25318,"none",10,8000),
  plot.pollen(dataset.25318,"m.avg",10,8000),
  plot.pollen(dataset.25318,"age.w",10,8000),
  plot.pollen(dataset.25318,"grim",10,8000),
  plot.pollen(dataset.25318,"shep",10,8000),
  ncol=5, nrow = 1, common.legend = T, legend = "right"
)

plot.pollen.25318



# comparison of result between diferent settings

dataset.25318.comparison.sample <- plot.comparison(dataset.25318,
                                            BIN = F, 
                                            BIN.size = 500, 
                                            Shiftbin = F, 
                                            N.shifts = 5, 
                                            rand = 1000, 
                                            interest.treshold =  8000)

dataset.25318.comparison.BIN <- plot.comparison(dataset.25318,
                                                BIN = T, 
                                                BIN.size = 500, 
                                                Shiftbin = F, 
                                                N.shifts = 5, 
                                                rand = 1000, 
                                                interest.treshold =  8000)

dataset.25318.comparison.BIN.shift <- plot.comparison(dataset.25318,
                                                BIN = T, 
                                                BIN.size = 500, 
                                                Shiftbin = T, 
                                                N.shifts = 5, 
                                                rand = 1000, 
                                                interest.treshold =  8000)

dataset.25318.comparison.sample
dataset.25318.comparison.BIN
dataset.25318.comparison.BIN.shift


# save.image("~/HOPE/GITHUB/RateOfChange/ENV_METHOD_20200319.RData")
# load("~/HOPE/GITHUB/RateOfChange/ENV_METHOD_20200319.RData")



# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())




RandomEnv(nforc = 1, time = 8e3:0, nprox=10)
test <- RandomProx()

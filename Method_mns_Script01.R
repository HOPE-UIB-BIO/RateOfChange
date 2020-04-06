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


# density values
low_diversity <- 5
high_diversity <- 50


plot.pollen <- function (data, sm.type, N.taxa, interest.treshold)
{
  Common.list <- data$filtered.counts %>%
    colSums() %>% 
    sort(decreasing = T) %>%
    .subset(.,1:N.taxa) %>%
    names() %>%
    sub("/",".",.) %>%
    sub("-",".",.) %>%
    sub(")",".",.) %>%
    sub(".\\(","..",.) 
    
  
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


plot.comparison <- function(data, BIN, BIN.size, Shiftbin, N.shifts, rand, standardise ,interest.treshold)
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


plot.time <- function(data, BIN=F, BIN.size=500, Shiftbin=F, N.shifts=5, rand=10)
{
  DF.performance <- data.frame(matrix(nrow = 20, ncol=5))
  names(DF.performance) <- c("smooth","DC","user","system","elapsed")
  DF.performance$smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4))
  DF.performance$DC <- c(rep(c("euc","euc.sd","chord","chisq"),5))
  
  for(i in 1:nrow(DF.performance))
  {
    a<- system.time(fc_ratepol( data.source.pollen =  data$filtered.counts,
                                data.source.age = data$list_ages,
                                sm.type = DF.performance$smooth[i],
                                N.points = 5,
                                range.age.max = 500, 
                                grim.N.max = 9,
                                BIN = BIN,
                                BIN.size = BIN.size,
                                Shiftbin = Shiftbin,
                                N.shifts = N.shifts,
                                rand = rand,
                                standardise = F,
                                DC = DF.performance$DC[i],
                                interest.treshold = 8000,
                                Debug = F))
    
    DF.performance$user[i] <- a[1]
    DF.performance$system[i] <- a[2]
    DF.performance$elapsed[i] <- a[3]
  }
  
  plot.fin <- DF.performance %>%
    ggplot(aes(y=elapsed, x=smooth, fill=DC))+
    geom_bar(stat="identity", position = "dodge", color="gray50")+
    ggtitle(paste(",N samples",nrow(data$filtered.counts),
                  ",N taxa",ncol(data$filtered.counts),
                  ",N.randomisation",rand,
                  ",BIN",BIN,
                  ",Shift",Shiftbin))+
    theme_classic()+
    coord_cartesian(ylim = c(0,60))+
    ylab("computation time (s)")
  return (plot.fin)
}


# -----------------------------------------
#
#       SIMULATION change 2000-3000 
#             low diversity
#
# -----------------------------------------

# ------------------------------
#             data
# ------------------------------

data_sim_ld_recent <-  fc_random_data(time = tibble_Europe2$list_ages[[2]]$ages$age,
                            nforc = 4, 
                            nprox = low_diversity,
                            manual.edit = T,
                            breaks=c(2000, 3000),
                            jitter = T,
                            rarity = T)

# ------------------------------
#          smoothing 
# ------------------------------

pollen_sim_ld_recent <- ggarrange(
  data_sim_ld_recent$list_ages$ages %>%
    filter(age <= 10000) %>%
    ggplot(aes(x=age))+
    geom_density(color="gray30",fill="gray80")+
    coord_flip(xlim = c(0,8000))+
    scale_x_continuous(trans = "reverse")+
    theme_classic()+xlab("Age")+ylab("Sample density")+
    ggtitle("Density of Samples") ,
  plot.pollen(data_sim_ld_recent,"none",low_diversity,8000),
  plot.pollen(data_sim_ld_recent,"m.avg",low_diversity,8000),
  plot.pollen(data_sim_ld_recent,"age.w",low_diversity,8000),
  plot.pollen(data_sim_ld_recent,"grim",low_diversity,8000),
  plot.pollen(data_sim_ld_recent,"shep",low_diversity,8000),
  ncol=6, nrow = 1, common.legend = T, legend = "none"
)

pollen_sim_ld_recent

ggsave("pollen_sim_ld_recent.pdf",
       plot = pollen_sim_ld_recent,
       width = 40, height = 15, units = "cm")

# ------------------------------
#   visual result comparison
# ------------------------------
visual_sim_ld_recent_sample <- plot.comparison(data_sim_ld_recent,
                                                   BIN = F, 
                                                   Shiftbin = F, 
                                                   rand = 1,
                                                   standardise = F,
                                                   interest.treshold =  8000)

ggsave("visual_sim_ld_recent_sample.pdf",
       plot = visual_sim_ld_recent_sample,
       width = 40, height = 25, units = "cm")

visual_sim_ld_recent_BIN <- plot.comparison(data_sim_ld_recent,
                                              BIN = T,
                                              BIN.size = 500,
                                              Shiftbin = F, 
                                              rand = 1, 
                                              standardise = F,
                                              interest.treshold =  8000)
ggsave("visual_sim_ld_recent_BIN.pdf",
       plot = visual_sim_ld_recent_BIN,
       width = 40, height = 25, units = "cm")


visual_sim_ld_recent_MW <- plot.comparison(data_sim_ld_recent,
                                                 BIN = T, 
                                                 BIN.size = 500,
                                                 Shiftbin = T, 
                                                 N.shifts = 5,
                                                 rand = 1, 
                                                 standardise = F,
                                                 interest.treshold =  8000)
ggsave("visual_sim_ld_recent_MW.pdf",
       plot = visual_sim_ld_recent_MW,
       width = 40, height = 25, units = "cm")

# ------------------------------
# statistical result comparison
# ------------------------------

perform_sim_ld_recent_sample <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                    nforc=4,
                                                    mean=100, 
                                                    sdev=.15, 
                                                    nprox=low_diversity, 
                                                    var=20,
                                                    range=15,
                                                    manual.edit = T,
                                                    breaks=c(2000,3000),
                                                    jitter = T,
                                                    rarity = T,
                                                    BIN=F,
                                                    Shiftbin=F,
                                                    rand.sets=100,
                                                    interest.treshold=8000)

ggsave("perform_sim_ld_recent_sample.pdf",
       plot = perform_sim_ld_recent_sample,
       width = 40, height = 25, units = "cm")



perform_sim_ld_recent_BIN <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                 nforc=4,
                                                 mean=100, 
                                                 sdev=.15, 
                                                 nprox=low_diversity, 
                                                 var=20,
                                                 range=15,
                                                 manual.edit = T,
                                                 breaks=c(2000,3000),
                                                 jitter = T,
                                                 rarity = T,
                                                 BIN=T,
                                                 BIN.size=500, 
                                                 Shiftbin=F,
                                                 rand.sets=100,
                                                 interest.treshold=8000)

ggsave("perform_sim_ld_recent_BIN.pdf",
       plot = perform_sim_ld_recent_BIN,
       width = 40, height = 25, units = "cm")



perform_sim_ld_recent_MW <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                           nforc=4,
                                                           mean=100, 
                                                           sdev=.15, 
                                                           nprox=10, 
                                                           var=20,
                                                           range=15,
                                                           manual.edit = T,
                                                           breaks=c(2000,3000),
                                                           jitter = T,
                                                           rarity = T,
                                                           BIN=T,
                                                           BIN.size=500, 
                                                           Shiftbin=T,
                                                           N.shifts=5,
                                                           rand.sets=100,
                                                           interest.treshold=8000)

ggsave("perform_sim_ld_recent_MW.pdf",
       plot = perform_sim_ld_recent_MW,
       width = 40, height = 25, units = "cm")


# -----------------------------------------
#
#       SIMULATION change 5500-6500 
#           low diversity
#
# -----------------------------------------

# ------------------------------
#               data
# ------------------------------
data_sim_ld_late <-  fc_random_data(time = tibble_Europe2$list_ages[[2]]$ages$age,
                               nforc = 4, 
                               nprox = low_diversity,
                               manual.edit = T,
                               breaks=c(5500, 6500),
                               jitter = T,
                               rarity = T)


# ------------------------------
#         smoothing
# ------------------------------

pollen_sim_ld_late <- ggarrange(
  data_sim_ld_late$list_ages$ages %>%
    filter(age <= 10000) %>%
    ggplot(aes(x=age))+
    geom_density(color="gray30",fill="gray80")+
    coord_flip(xlim = c(0,8000))+
    scale_x_continuous(trans = "reverse")+
    theme_classic()+xlab("Age")+ylab("Sample density")+
    ggtitle("Density of Samples") ,
  plot.pollen(data_sim_ld_late,"none",low_diversity,8000),
  plot.pollen(data_sim_ld_late,"m.avg",low_diversity,8000),
  plot.pollen(data_sim_ld_late,"age.w",low_diversity,8000),
  plot.pollen(data_sim_ld_late,"grim",low_diversity,8000),
  plot.pollen(data_sim_ld_late,"shep",low_diversity,8000),
  ncol=6, nrow = 1, common.legend = T, legend = "none"
)

pollen_sim_ld_late


ggsave("pollen_sim_ld_late.pdf",
       plot = pollen_sim_ld_late,
       width = 40, height = 15, units = "cm")


# ------------------------------
#   visual result comparison
# ------------------------------


visual_sim_ld_late_sample <- plot.comparison(data_sim_ld_late,
                                                    BIN = F, 
                                                    Shiftbin = F, 
                                                    rand = 1, 
                                                    interest.treshold =  8000)

ggsave("visual_sim_ld_late_sample.pdf",
       plot = visual_sim_ld_late_sample,
       width = 40, height = 25, units = "cm")

visual_sim_ld_late_BIN <- plot.comparison(data_sim_ld_late,
                                                 BIN = T,
                                                 BIN.size = 500,
                                                 Shiftbin = F, 
                                                 rand = 1, 
                                                 interest.treshold =  8000)
ggsave("visual_sim_ld_late_BIN.pdf",
       plot = visual_sim_ld_late_BIN,
       width = 40, height = 25, units = "cm")


visual_sim_ld_late_MW <- plot.comparison(data_sim_ld_late,
                                                       BIN = T, 
                                                       BIN.size = 500,
                                                       Shiftbin = T, 
                                                       N.shifts = 5,
                                                       rand = 1, 
                                                       interest.treshold =  8000)
ggsave("visual_sim_ld_late_MW.pdf",
       plot = visual_sim_ld_late_MW,
       width = 40, height = 25, units = "cm")


# ------------------------------
# statistical result comparison
# ------------------------------

perform_sim_ld_late_sample <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                       nforc=4,
                                                       mean=100, 
                                                       sdev=.15, 
                                                       nprox=low_diversity, 
                                                       var=20,
                                                       range=15,
                                                       manual.edit = T,
                                                       breaks=c(5500,6500),
                                                       jitter = T,
                                                       rarity = T,
                                                       BIN=F,
                                                       Shiftbin=F,
                                                       rand.sets=10,
                                                       interest.treshold=8000)

ggsave("perform_sim_ld_late_sample.pdf",
       plot = perform_sim_ld_late_sample,
       width = 40, height = 25, units = "cm")



perform_sim_ld_late_BIN <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                    nforc=4,
                                                    mean=100, 
                                                    sdev=.15, 
                                                    nprox=low_diversity, 
                                                    var=20,
                                                    range=15,
                                                    manual.edit = T,
                                                    breaks=c(5500,6500),
                                                    jitter = T,
                                                    rarity = T,
                                                    BIN=T,
                                                    BIN.size=500, 
                                                    Shiftbin=F,
                                                    rand.sets=10,
                                                    interest.treshold=8000)

ggsave("perform_sim_ld_late_BIN.pdf",
       plot = perform_sim_ld_late_BIN,
       width = 40, height = 25, units = "cm")



perform_sim_ld_late_MW <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                           nforc=4,
                                                           mean=100, 
                                                           sdev=.15, 
                                                           nprox=low_diversity, 
                                                           var=20,
                                                           range=15,
                                                           manual.edit = T,
                                                           breaks=c(5500,6500),
                                                           jitter = T,
                                                           rarity = T,
                                                           BIN=T,
                                                           BIN.size=500, 
                                                           Shiftbin=T,
                                                           N.shifts=5,
                                                           rand.sets=10,
                                                           interest.treshold=8000)

ggsave("perform_sim_ld_late_MW.pdf",
       plot = perform_sim_ld_late_MW,
       width = 40, height = 25, units = "cm")



# -----------------------------------------
#
#       SIMULATION change 2000-3000 
#             high diversity
#
# -----------------------------------------

# ------------------------------
#             data
# ------------------------------

data_sim_hd_recent <-  fc_random_data(time = tibble_Europe2$list_ages[[2]]$ages$age,
                               nforc = 4, 
                               nprox = high_diversity,
                               manual.edit = T,
                               breaks=c(2000, 3000),
                               jitter = T,
                               rarity = T)


# ------------------------------
#          smoothing 
# ------------------------------

pollen_sim_hd_recent <- ggarrange(
  data_sim_hd_recent$list_ages$ages %>%
    filter(age <= 10000) %>%
    ggplot(aes(x=age))+
    geom_density(color="gray30",fill="gray80")+
    coord_flip(xlim = c(0,8000))+
    scale_x_continuous(trans = "reverse")+
    theme_classic()+xlab("Age")+ylab("Sample density")+
    ggtitle("Density of Samples") ,
  plot.pollen(data_sim_hd_recent,"none",high_diversity,8000),
  plot.pollen(data_sim_hd_recent,"m.avg",high_diversity,8000),
  plot.pollen(data_sim_hd_recent,"age.w",high_diversity,8000),
  plot.pollen(data_sim_hd_recent,"grim",high_diversity,8000),
  plot.pollen(data_sim_hd_recent,"shep",high_diversity,8000),
  ncol=6, nrow = 1, common.legend = T, legend = "none"
)

pollen_sim_hd_recent

ggsave("pollen_sim_hd_recent.pdf",
       plot = pollen_sim_hd_recent,
       width = 40, height = 15, units = "cm")

# ------------------------------
#   visual result comparison
# ------------------------------
visual_sim_hd_recent_sample <- plot.comparison(data_sim_hd_recent,
                                                    BIN = F, 
                                                    Shiftbin = F, 
                                                    rand = 1,
                                                    standardise = F,
                                                    interest.treshold =  8000)

ggsave("visual_sim_hd_recent_sample.pdf",
       plot = visual_sim_hd_recent_sample,
       width = 40, height = 25, units = "cm")

visual_sim_hd_recent_BIN <- plot.comparison(data_sim_hd_recent,
                                                 BIN = T,
                                                 BIN.size = 500,
                                                 Shiftbin = F, 
                                                 rand = 1, 
                                                 standardise = F,
                                                 interest.treshold =  8000)
ggsave("visual_sim_hd_recent_BIN.pdf",
       plot = visual_sim_hd_recent_BIN,
       width = 40, height = 25, units = "cm")


visual_sim_hd_recent_MW <- plot.comparison(data_sim_hd_recent,
                                                       BIN = T, 
                                                       BIN.size = 500,
                                                       Shiftbin = T, 
                                                       N.shifts = 5,
                                                       rand = 1, 
                                                       standardise = F,
                                                       interest.treshold =  8000)
ggsave("visual_sim_hd_recent_MW.pdf",
       plot = visual_sim_hd_recent_MW,
       width = 40, height = 25, units = "cm")

# ------------------------------
# statistical result comparison
# ------------------------------

perform_sim_hd_recent_sample <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                       nforc=4,
                                                       mean=100, 
                                                       sdev=.15, 
                                                       nprox=high_diversity, 
                                                       var=20,
                                                       range=15,
                                                       manual.edit = T,
                                                       breaks=c(2000,3000),
                                                       jitter = T,
                                                       rarity = T,
                                                       BIN=F,
                                                       Shiftbin=F,
                                                       rand.sets=100,
                                                       interest.treshold=8000)

ggsave("perform_sim_hd_recent_sample.pdf",
       plot = perform_sim_hd_recent_sample,
       width = 40, height = 25, units = "cm")



perform_sim_hd_recent_BIN <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                    nforc=4,
                                                    mean=100, 
                                                    sdev=.15, 
                                                    nprox=high_diversity, 
                                                    var=20,
                                                    range=15,
                                                    manual.edit = T,
                                                    breaks=c(2000,3000),
                                                    jitter = T,
                                                    rarity = T,
                                                    BIN=T,
                                                    BIN.size=500, 
                                                    Shiftbin=F,
                                                    rand.sets=100,
                                                    interest.treshold=8000)

ggsave("perform_sim_hd_recent_BIN.pdf",
       plot = perform_sim_hd_recent_BIN,
       width = 40, height = 25, units = "cm")



perform_sim_hd_recent_MW <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                           nforc=4,
                                                           mean=100, 
                                                           sdev=.15, 
                                                           nprox=high_diversity, 
                                                           var=20,
                                                           range=15,
                                                           manual.edit = T,
                                                           breaks=c(2000,3000),
                                                           jitter = T,
                                                           rarity=T,
                                                           BIN=T,
                                                           BIN.size=500, 
                                                           Shiftbin=T,
                                                           N.shifts=5,
                                                           rand.sets=10,
                                                           interest.treshold=8000)

ggsave("perform_sim_hd_recent_MW.pdf",
       plot = perform_sim_hd_recent_MW,
       width = 40, height = 25, units = "cm")



# -----------------------------------------
#
#       SIMULATION change 5500-6500 
#           high diversity
#
# -----------------------------------------

# ------------------------------
#               data
# ------------------------------
data_sim_hd_late <-  fc_random_data(time = tibble_Europe2$list_ages[[2]]$ages$age,
                                    nforc = 4, 
                                    nprox = high_diversity,
                                    manual.edit = T,
                                    breaks=c(5500, 6500),
                                    jitter = T)


# ------------------------------
#         smoothing
# ------------------------------

pollen_sim_hd_late <- ggarrange(
  data_sim_hd_late$list_ages$ages %>%
    filter(age <= 10000) %>%
    ggplot(aes(x=age))+
    geom_density(color="gray30",fill="gray80")+
    coord_flip(xlim = c(0,8000))+
    scale_x_continuous(trans = "reverse")+
    theme_classic()+xlab("Age")+ylab("Sample density")+
    ggtitle("Density of Samples") ,
  plot.pollen(data_sim_hd_late,"none",high_diversity,8000),
  plot.pollen(data_sim_hd_late,"m.avg",high_diversity,8000),
  plot.pollen(data_sim_hd_late,"age.w",high_diversity,8000),
  plot.pollen(data_sim_hd_late,"grim",high_diversity,8000),
  plot.pollen(data_sim_hd_late,"shep",high_diversity,8000),
  ncol=6, nrow = 1, common.legend = T, legend = "none"
)

pollen_sim_hd_late


ggsave("pollen_sim_hd_late.pdf",
       plot = pollen_sim_hd_late,
       width = 40, height = 15, units = "cm")


# ------------------------------
#   visual result comparison
# ------------------------------


visual_sim_hd_late_sample <- plot.comparison(data_sim_hd_late,
                                             BIN = F, 
                                             Shiftbin = F, 
                                             rand = 1, 
                                             interest.treshold =  8000)

ggsave("visual_sim_hd_late_sample.pdf",
       plot = visual_sim_hd_late_sample,
       width = 40, height = 25, units = "cm")

visual_sim_hd_late_BIN <- plot.comparison(data_sim_hd_late,
                                          BIN = T,
                                          BIN.size = 500,
                                          Shiftbin = F, 
                                          rand = 1, 
                                          interest.treshold =  8000)
ggsave("visual_sim_hd_late_BIN.pdf",
       plot = visual_sim_hd_late_BIN,
       width = 40, height = 25, units = "cm")


visual_sim_hd_late_MW <- plot.comparison(data_sim_hd_late,
                                         BIN = T, 
                                         BIN.size = 500,
                                         Shiftbin = T, 
                                         N.shifts = 5,
                                         rand = 1, 
                                         interest.treshold =  8000)
ggsave("visual_sim_hd_late_MW.pdf",
       plot = visual_sim_hd_late_MW,
       width = 40, height = 25, units = "cm")


# ------------------------------
# statistical result comparison
# ------------------------------

perform_sim_hd_late_sample <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                                  nforc=4,
                                                  mean=100, 
                                                  sdev=.15, 
                                                  nprox=high_diversity, 
                                                  var=20,
                                                  range=15,
                                                  manual.edit = T,
                                                  breaks=c(5500,6500),
                                                  jitter = T,
                                                  rarity = T,
                                                  BIN=F,
                                                  Shiftbin=F,
                                                  rand.sets=10,
                                                  interest.treshold=8000)

ggsave("perform_sim_hd_late_sample.pdf",
       plot = perform_sim_hd_late_sample,
       width = 40, height = 25, units = "cm")



perform_sim_hd_late_BIN <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                               nforc=4,
                                               mean=100, 
                                               sdev=.15, 
                                               nprox=high_diversity, 
                                               var=20,
                                               range=15,
                                               manual.edit = T,
                                               breaks=c(5500,6500),
                                               jitter = T,
                                               rarity = T,
                                               BIN=T,
                                               BIN.size=500, 
                                               Shiftbin=F,
                                               rand.sets=10,
                                               interest.treshold=8000)

ggsave("perform_sim_hd_late_BIN.pdf",
       plot = perform_sim_hd_late_BIN,
       width = 40, height = 25, units = "cm")



perform_sim_hd_late_MW <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
                                              nforc=4,
                                              mean=100, 
                                              sdev=.15, 
                                              nprox=high_diversity, 
                                              var=20,
                                              range=15,
                                              manual.edit = T,
                                              breaks=c(5500,6500),
                                              jitter = T,
                                              rarity = T,
                                              BIN=T,
                                              BIN.size=500, 
                                              Shiftbin=T,
                                              N.shifts=5,
                                              rand.sets=10,
                                              interest.treshold=8000)

ggsave("perform_sim_hd_late_MW.pdf",
       plot = perform_sim_hd_late_MW,
       width = 40, height = 25, units = "cm")


# -----------------------------------------
#
#                 17334
#
# -----------------------------------------


data_17334 <- list(dataset.id = tibble_Europe2$dataset.id[[2]],
                      filtered.counts = tibble_Europe2$filtered.counts[[2]],
                      list_ages = tibble_Europe2$list_ages[[2]])
  
# POllen graph

pollen_17334 <- ggarrange(
  data_17334$list_ages$ages %>%
    filter(age <= 10000) %>%
    ggplot(aes(x=age))+
    geom_density(color="gray30",fill="gray80")+
    coord_flip(xlim = c(0,8000))+
    scale_x_continuous(trans = "reverse")+
    theme_classic()+xlab("Age")+ylab("Sample density")+
    ggtitle("Density of Samples"),
  plot.pollen(data_17334,"none",10,8000),
  plot.pollen(data_17334,"m.avg",10,8000),
  plot.pollen(data_17334,"age.w",10,8000),
  plot.pollen(data_17334,"grim",10,8000),
  plot.pollen(data_17334,"shep",10,8000),
  ncol=6, nrow = 1, common.legend = T, legend = "right"
)

pollen_17334

ggsave("pollen_17334.pdf",
       plot = pollen_17334,
       width = 40, height = 15, units = "cm")


# ------------------------------
#   visual result comparison
# ------------------------------
visual_17334_sample <- plot.comparison(data_17334,
                                        BIN = F, 
                                        Shiftbin = F, 
                                        rand = 1000,
                                        standardise = T,
                                        interest.treshold =  8000)

ggsave("visual_17334_sample.pdf",
       plot = visual_17334_sample,
       width = 40, height = 25, units = "cm")

visual_17334_BIN <- plot.comparison(data_17334,
                                                 BIN = T,
                                                 BIN.size = 500,
                                                 Shiftbin = F, 
                                                 rand = 1000,
                                                 standardise = T,
                                                 interest.treshold =  8000)
ggsave("visual_17334_BIN",
       plot = visual_17334_BIN,
       width = 40, height = 25, units = "cm")


visual_17334_MW <- plot.comparison(data_17334,
                                                       BIN = T, 
                                                       BIN.size = 500,
                                                       Shiftbin = T, 
                                                       N.shifts = 5,
                                                       rand = 1000,
                                                       standardise = T,
                                                       interest.treshold =  8000)
ggsave("visual_17334_MW.pdf",
       plot = visual_17334_MW,
       width = 40, height = 25, units = "cm")

# ------------------------------
#   computation time compariosn
# ------------------------------

time_17334_sample <- plot.time(data_17334, BIN=F, Shiftbin = F, rand = 1000)

time_17334_BIN <- plot.time(data_17334, BIN=T,BIN.size = 500, Shiftbin = F, rand = 1000)

time_17334_MW <- plot.time(data_17334, BIN=T,BIN.size = 500, Shiftbin = T, N.shifts = 5, rand = 1000)

time_17334_sum <- ggarrange(
  time_17334_sample,
  time_17334_BIN,
  time_17334_MW,
  nrow = 3, ncol = 1, common.legend = T, legend = "right"
)

ggsave("time_17334_sum.pdf",
       plot = time_17334_sum,
       height = 25, width = 25, units="cm")


# -----------------------------------------

# save.image("~/HOPE/GITHUB/RateOfChange/ENV_METHOD_20200319.RData")
# load("~/HOPE/GITHUB/RateOfChange/ENV_METHOD_20200319.RData")



# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())



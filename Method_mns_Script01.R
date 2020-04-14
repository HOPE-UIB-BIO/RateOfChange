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
# download.file("https://www.dropbox.com/s/3hp7rv03mkg4pjz/tibble_Europe_filtered05.03.20.RData?dl=1","~/input/DATA/tibble_Europe_filtered05.03.20.RData")

setwd("~/GITHUB/RateOfChange")

load("~/DATA/input/tibble_Europe_filtered05.03.20.RData")

files.sources <- list.files("~/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)


glimpse(tibble_Europe2)



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
  
  plot.data <- data.ext$Pollen %>%
    select(Common.list) %>%
    rownames_to_column() %>%
    pivot_longer(., cols = c(Common.list)) %>%
    rename(sample.id = rowname) %>%
    inner_join(.,data.ext$Age, by="sample.id")
  
  plot.p <-  plot.data %>%
    ggplot(aes( y=value, 
                x= age))+
    theme_classic()+
    scale_x_continuous(trans = "reverse")+
    geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
                color="gray20", alpha=1/5, size=0.1)+
    xlab("Age")+ylab("Pollen (%)")+
    coord_flip(xlim=c(interest.treshold,0), ylim = c(0,1))+
    ggtitle(paste("smoothing",sm.type))  
  return (list(data=plot.data,plot=plot.p))
}

plot.smooting <- function(data)
{
res.temp<-  ggarrange(
    data$list_ages$ages %>%
      filter(age <= age_lim) %>%
      ggplot(aes(x=age))+
      geom_density(color="gray30",fill="gray80")+
      scale_x_reverse()+
      coord_flip(xlim = c(age_lim,0))+
      theme_classic()+xlab("Age")+ylab("Sample density")+
      ggtitle("Density of Samples") ,
    plot.pollen(data,"none",low_diversity,age_lim)$plot,
    plot.pollen(data,"m.avg",low_diversity,age_lim)$plot,
    plot.pollen(data,"age.w",low_diversity,age_lim)$plot,
    plot.pollen(data,"grim",low_diversity,age_lim)$plot,
    plot.pollen(data,"shep",low_diversity,age_lim)$plot,
    ncol=6, nrow = 1, common.legend = T, legend = "none")  
return(res.temp)
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
    
    data.temp.sum<- data.frame(data.temp,performance.smooth[i],performance.DC[i])
    names(data.temp.sum) <- c(names(data.temp),"smooth","DC")
    
    if(i == 1){
      perfomance.tibble <- data.temp.sum} else {
        perfomance.tibble <- rbind(perfomance.tibble,data.temp.sum)
      }
    
    
    performance.list.plot[[i]]<-ggplot(data=data.temp, 
                                       aes(y=RUN.RoC, 
                                           x= RUN.Age.Pos))+
      theme_classic()+
      scale_x_continuous(trans = "reverse")+
      coord_flip(xlim=c(age.max,0), ylim = c(0,roc.max))+
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
  
  return(list(data=perfomance.tibble, plot=p.all))
}


plot.time <- function(data, BIN=F, BIN.size=500, Shiftbin=F, N.shifts=5, rand=100)
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
                                interest.treshold = age_lim,
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
    coord_cartesian(ylim = c(0,120))+
    ylab("computation time (s)")
  return (list(data=DF.performance, plot=plot.fin))
}

# ----------------------------------------------
#             DEFINITION OF VARIBLES
# ----------------------------------------------

# Number of simulated enviromental variables
N_env <- 4

# diversity of pollen taxat in simulated data
low_diversity <- 5
high_diversity <- 50

# position of the enviromental change in the sequence 
breaks_recent <- c(2000, 3000)
breaks_late <- c(5500, 6500)

# Number of repetion in simulation pollen data
N_rep <- 100

# template of time sequence with uneven distribution of points
time_seq <- tibble_Europe2$list_ages[[2]]$ages$age

# maximal time of focus 
age_lim <- max(time_seq)

# -----------------------------------------
#
#           SIMULATION DATA 
#
# -----------------------------------------

# low diversity - recent
data_sim_ld_recent <-  fc_random_data(time = time_seq,
                                      nforc = N_env, 
                                      nprox = low_diversity,
                                      manual.edit = T,
                                      breaks=breaks_recent,
                                      jitter = T,
                                      rarity = T)


pollen_sim_ld_recent <- plot.smooting(data_sim_ld_recent)

pollen_sim_ld_recent

ggsave("~RESULTS/Methods/pollen_sim_ld_recent.pdf",
       plot = pollen_sim_ld_recent,
       width = 40, height = 15, units = "cm")

# low diversity - late

data_sim_ld_late <-  fc_random_data(time = time_seq,
                                    nforc = N_env, 
                                    nprox = low_diversity,
                                    manual.edit = T,
                                    breaks=breaks_late,
                                    jitter = T,
                                    rarity = T)


pollen_sim_ld_late <- plot.smooting(data_sim_ld_late)

pollen_sim_ld_late


ggsave("~RESULTS/Methods/pollen_sim_ld_late.pdf",
       plot = pollen_sim_ld_late,
       width = 40, height = 15, units = "cm")


# hight diversity - recent

data_sim_hd_recent <-  fc_random_data(time = time_seq,
                                      nforc = N_env, 
                                      nprox = high_diversity,
                                      manual.edit = T,
                                      breaks=breaks_recent,
                                      jitter = T,
                                      rarity = T)

pollen_sim_hd_recent <- plot.smooting(data_sim_hd_recent)

pollen_sim_hd_recent

ggsave("~RESULTS/Methods/pollen_sim_hd_recent.pdf",
       plot = pollen_sim_hd_recent,
       width = 40, height = 15, units = "cm")


# high diversity - late

data_sim_hd_late <-  fc_random_data(time = time_seq,
                                    nforc = N_env, 
                                    nprox = high_diversity,
                                    manual.edit = T,
                                    breaks=breaks_late,
                                    jitter = T)


pollen_sim_hd_late <- plot.smooting(data_sim_hd_late)

pollen_sim_hd_late

ggsave("~RESULTS/Methods/pollen_sim_hd_late.pdf",
       plot = pollen_sim_hd_late,
       width = 40, height = 15, units = "cm")


# -----------------------------------------
#
#             VISUAL RESULTS
# 
# -----------------------------------------

# low diversity - recent

visual_sim_ld_recent_sample <- plot.comparison(data_sim_ld_recent,
                                               BIN = F, 
                                               Shiftbin = F, 
                                               rand = 1,
                                               standardise = F,
                                               interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_ld_recent_sample.pdf",
       plot = visual_sim_ld_recent_sample$plot,
       width = 40, height = 25, units = "cm")

visual_sim_ld_recent_BIN <- plot.comparison(data_sim_ld_recent,
                                            BIN = T,
                                            BIN.size = 500,
                                            Shiftbin = F, 
                                            rand = 1, 
                                            standardise = F,
                                            interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_ld_recent_BIN.pdf",
       plot = visual_sim_ld_recent_BIN$plot,
       width = 40, height = 25, units = "cm")


visual_sim_ld_recent_MW <- plot.comparison(data_sim_ld_recent,
                                           BIN = T, 
                                           BIN.size = 500,
                                           Shiftbin = T, 
                                           N.shifts = 5,
                                           rand = 1, 
                                           standardise = F,
                                           interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_ld_recent_MW.pdf",
       plot = visual_sim_ld_recent_MW$plot,
       width = 40, height = 25, units = "cm")



# low diversity late


visual_sim_ld_late_sample <- plot.comparison(data_sim_ld_late,
                                             BIN = F, 
                                             Shiftbin = F, 
                                             rand = 1, 
                                             interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_ld_late_sample.pdf",
       plot = visual_sim_ld_late_sample$plot,
       width = 40, height = 25, units = "cm")

visual_sim_ld_late_BIN <- plot.comparison(data_sim_ld_late,
                                          BIN = T,
                                          BIN.size = 500,
                                          Shiftbin = F, 
                                          rand = 1, 
                                          interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_ld_late_BIN.pdf",
       plot = visual_sim_ld_late_BIN$plot,
       width = 40, height = 25, units = "cm")


visual_sim_ld_late_MW <- plot.comparison(data_sim_ld_late,
                                         BIN = T, 
                                         BIN.size = 500,
                                         Shiftbin = T, 
                                         N.shifts = 5,
                                         rand = 1, 
                                         interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_ld_late_MW.pdf",
       plot = visual_sim_ld_late_MW$plot,
       width = 40, height = 25, units = "cm")

# high diveristy - recent
visual_sim_hd_recent_sample <- plot.comparison(data_sim_hd_recent,
                                               BIN = F, 
                                               Shiftbin = F, 
                                               rand = 1,
                                               standardise = F,
                                               interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_hd_recent_sample.pdf",
       plot = visual_sim_hd_recent_sample$plot,
       width = 40, height = 25, units = "cm")

visual_sim_hd_recent_BIN <- plot.comparison(data_sim_hd_recent,
                                            BIN = T,
                                            BIN.size = 500,
                                            Shiftbin = F, 
                                            rand = 1, 
                                            standardise = F,
                                            interest.treshold =  age_lim)
ggsave("~RESULTS/Methods/visual_sim_hd_recent_BIN.pdf",
       plot = visual_sim_hd_recent_BIN$plot,
       width = 40, height = 25, units = "cm")


visual_sim_hd_recent_MW <- plot.comparison(data_sim_hd_recent,
                                           BIN = T, 
                                           BIN.size = 500,
                                           Shiftbin = T, 
                                           N.shifts = 5,
                                           rand = 1, 
                                           standardise = F,
                                           interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_hd_recent_MW.pdf",
       plot = visual_sim_hd_recent_MW$plot,
       width = 40, height = 25, units = "cm")


# hight diversity - late


visual_sim_hd_late_sample <- plot.comparison(data_sim_hd_late,
                                             BIN = F, 
                                             Shiftbin = F, 
                                             rand = 1, 
                                             interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_hd_late_sample.pdf",
       plot = visual_sim_hd_late_sample$plot,
       width = 40, height = 25, units = "cm")

visual_sim_hd_late_BIN <- plot.comparison(data_sim_hd_late,
                                          BIN = T,
                                          BIN.size = 500,
                                          Shiftbin = F, 
                                          rand = 1, 
                                          interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_hd_late_BIN.pdf",
       plot = visual_sim_hd_late_BIN$plot,
       width = 40, height = 25, units = "cm")


visual_sim_hd_late_MW <- plot.comparison(data_sim_hd_late,
                                         BIN = T, 
                                         BIN.size = 500,
                                         Shiftbin = T, 
                                         N.shifts = 5,
                                         rand = 1, 
                                         interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_sim_hd_late_MW.pdf",
       plot = visual_sim_hd_late_MW$plot,
       width = 40, height = 25, units = "cm")



# -----------------------------------------
#
#             STATISTICAL COMP
# 
# -----------------------------------------

# low diversity - recent

perform_sim_ld_recent_sample <- fc_random_data_test(time= time_seq,
                                                    nforc=4,
                                                    mean=100, 
                                                    sdev=.15, 
                                                    nprox=low_diversity, 
                                                    var=20,
                                                    range=15,
                                                    manual.edit = T,
                                                    breaks=breaks_recent,
                                                    jitter = T,
                                                    rarity = T,
                                                    BIN=F,
                                                    Shiftbin=F,
                                                    rand.sets=N_rep,
                                                    interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_ld_recent_sample.pdf",
       plot = perform_sim_ld_recent_sample$plot,
       width = 40, height = 25, units = "cm")



perform_sim_ld_recent_BIN <- fc_random_data_test(time= time_seq,
                                                 nforc=N_env,
                                                 mean=100, 
                                                 sdev=.15, 
                                                 nprox=low_diversity, 
                                                 var=20,
                                                 range=15,
                                                 manual.edit = T,
                                                 breaks=breaks_recent,
                                                 jitter = T,
                                                 rarity = T,
                                                 BIN=T,
                                                 BIN.size=500, 
                                                 Shiftbin=F,
                                                 rand.sets=N_rep,
                                                 interest.treshold=age_lim)




ggsave("~RESULTS/Methods/perform_sim_ld_recent_BIN.pdf",
       plot = perform_sim_ld_recent_BIN$plot,
       width = 40, height = 25, units = "cm")



perform_sim_ld_recent_MW <- fc_random_data_test(time= time_seq,
                                                nforc=N_env,
                                                mean=100, 
                                                sdev=.15, 
                                                nprox=10, 
                                                var=20,
                                                range=15,
                                                manual.edit = T,
                                                breaks=breaks_recent,
                                                jitter = T,
                                                rarity = T,
                                                BIN=T,
                                                BIN.size=500, 
                                                Shiftbin=T,
                                                N.shifts=5,
                                                rand.sets=N_rep,
                                                interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_ld_recent_MW.pdf",
       plot = perform_sim_ld_recent_MW$plot,
       width = 40, height = 25, units = "cm")



# low diversity - late 

perform_sim_ld_late_sample <- fc_random_data_test(time= time_seq,
                                                  nforc=N_env,
                                                  mean=100, 
                                                  sdev=.15, 
                                                  nprox=low_diversity, 
                                                  var=20,
                                                  range=15,
                                                  manual.edit = T,
                                                  breaks=breaks_late,
                                                  jitter = T,
                                                  rarity = T,
                                                  BIN=F,
                                                  Shiftbin=F,
                                                  rand.sets=N_rep,
                                                  interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_ld_late_sample.pdf",
       plot = perform_sim_ld_late_sample$plot,
       width = 40, height = 25, units = "cm")



perform_sim_ld_late_BIN <- fc_random_data_test(time= time_seq,
                                               nforc=N_env,
                                               mean=100, 
                                               sdev=.15, 
                                               nprox=low_diversity, 
                                               var=20,
                                               range=15,
                                               manual.edit = T,
                                               breaks=breaks_late,
                                               jitter = T,
                                               rarity = T,
                                               BIN=T,
                                               BIN.size=500, 
                                               Shiftbin=F,
                                               rand.sets=N_rep,
                                               interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_ld_late_BIN.pdf",
       plot = perform_sim_ld_late_BIN$plot,
       width = 40, height = 25, units = "cm")



perform_sim_ld_late_MW <- fc_random_data_test(time= time_seq,
                                              nforc=N_env,
                                              mean=100, 
                                              sdev=.15, 
                                              nprox=low_diversity, 
                                              var=20,
                                              range=15,
                                              manual.edit = T,
                                              breaks=breaks_late,
                                              jitter = T,
                                              rarity = T,
                                              BIN=T,
                                              BIN.size=500, 
                                              Shiftbin=T,
                                              N.shifts=5,
                                              rand.sets=N_rep,
                                              interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_ld_late_MW.pdf",
       plot = perform_sim_ld_late_MW$plot,
       width = 40, height = 25, units = "cm")



# high diversity  - recent

perform_sim_hd_recent_sample <- fc_random_data_test(time= time_seq,
                                                    nforc=N_env,
                                                    mean=100, 
                                                    sdev=.15, 
                                                    nprox=high_diversity, 
                                                    var=20,
                                                    range=15,
                                                    manual.edit = T,
                                                    breaks=breaks_recent,
                                                    jitter = T,
                                                    rarity = T,
                                                    BIN=F,
                                                    Shiftbin=F,
                                                    rand.sets=N_rep,
                                                    interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_hd_recent_sample.pdf",
       plot = perform_sim_hd_recent_sample$plot,
       width = 40, height = 25, units = "cm")



perform_sim_hd_recent_BIN <- fc_random_data_test(time= time_seq,
                                                 nforc=N_env,
                                                 mean=100, 
                                                 sdev=.15, 
                                                 nprox=high_diversity, 
                                                 var=20,
                                                 range=15,
                                                 manual.edit = T,
                                                 breaks=breaks_recent,
                                                 jitter = T,
                                                 rarity = T,
                                                 BIN=T,
                                                 BIN.size=500, 
                                                 Shiftbin=F,
                                                 rand.sets=N_rep,
                                                 interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_hd_recent_BIN.pdf",
       plot = perform_sim_hd_recent_BIN$plot,
       width = 40, height = 25, units = "cm")



perform_sim_hd_recent_MW <- fc_random_data_test(time= time_seq,
                                                nforc=N_env,
                                                mean=100, 
                                                sdev=.15, 
                                                nprox=high_diversity, 
                                                var=20,
                                                range=15,
                                                manual.edit = T,
                                                breaks=breaks_recent,
                                                jitter = T,
                                                rarity=T,
                                                BIN=T,
                                                BIN.size=500, 
                                                Shiftbin=T,
                                                N.shifts=5,
                                                rand.sets=N_rep,
                                                interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_hd_recent_MW.pdf",
       plot = perform_sim_hd_recent_MW$plot,
       width = 40, height = 25, units = "cm")



# high diversity  - late

perform_sim_hd_late_sample <- fc_random_data_test(time= time_seq,
                                                  nforc=N_env,
                                                  mean=100, 
                                                  sdev=.15, 
                                                  nprox=high_diversity, 
                                                  var=20,
                                                  range=15,
                                                  manual.edit = T,
                                                  breaks=breaks_late,
                                                  jitter = T,
                                                  rarity = T,
                                                  BIN=F,
                                                  Shiftbin=F,
                                                  rand.sets=N_rep,
                                                  interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_hd_late_sample.pdf",
       plot = perform_sim_hd_late_sample$plot,
       width = 40, height = 25, units = "cm")



perform_sim_hd_late_BIN <- fc_random_data_test(time= time_seq,
                                               nforc=N_env,
                                               mean=100, 
                                               sdev=.15, 
                                               nprox=high_diversity, 
                                               var=20,
                                               range=15,
                                               manual.edit = T,
                                               breaks=breaks_late,
                                               jitter = T,
                                               rarity = T,
                                               BIN=T,
                                               BIN.size=500, 
                                               Shiftbin=F,
                                               rand.sets=10,
                                               interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_hd_late_BIN.pdf",
       plot = perform_sim_hd_late_BIN$plot,
       width = 40, height = 25, units = "cm")



perform_sim_hd_late_MW <- fc_random_data_test(time= time_seq,
                                              nforc=N_env,
                                              mean=100, 
                                              sdev=.15, 
                                              nprox=high_diversity, 
                                              var=20,
                                              range=15,
                                              manual.edit = T,
                                              breaks=breaks_late,
                                              jitter = T,
                                              rarity = T,
                                              BIN=T,
                                              BIN.size=500, 
                                              Shiftbin=T,
                                              N.shifts=5,
                                              rand.sets=N_rep,
                                              interest.treshold=age_lim)

ggsave("~RESULTS/Methods/perform_sim_hd_late_MW.pdf",
       plot = perform_sim_hd_late_MW$plot,
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
    filter(age <= age_lim) %>%
    ggplot(aes(x=age))+
    geom_density(color="gray30",fill="gray80")+
    coord_flip(xlim = c(age_lim,0))+
    scale_x_continuous(trans = "reverse")+
    theme_classic()+xlab("Age")+ylab("Sample density")+
    ggtitle("Density of Samples"),
  plot.pollen(data_17334,"none",10,age_lim)$plot,
  plot.pollen(data_17334,"m.avg",10,age_lim)$plot,
  plot.pollen(data_17334,"age.w",10,age_lim)$plot,
  plot.pollen(data_17334,"grim",10,age_lim)$plot,
  plot.pollen(data_17334,"shep",10,age_lim)$plot,
  ncol=6, nrow = 1, common.legend = T, legend = "right"
)

pollen_17334

ggsave("~RESULTS/Methods/pollen_17334.pdf",
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
                                        interest.treshold =  age_lim)

ggsave("~RESULTS/Methods/visual_17334_sample.pdf",
       plot = visual_17334_sample$plot,
       width = 40, height = 25, units = "cm")

visual_17334_BIN <- plot.comparison(data_17334,
                                                 BIN = T,
                                                 BIN.size = 500,
                                                 Shiftbin = F, 
                                                 rand = 1000,
                                                 standardise = T,
                                                 interest.treshold =  age_lim)
ggsave("~RESULTS/Methods/visual_17334_BIN.pdf",
       plot = visual_17334_BIN$plot,
       width = 40, height = 25, units = "cm")


visual_17334_MW <- plot.comparison(data_17334,
                                                       BIN = T, 
                                                       BIN.size = 500,
                                                       Shiftbin = T, 
                                                       N.shifts = 5,
                                                       rand = 1000,
                                                       standardise = T,
                                                       interest.treshold =  age_lim)
ggsave("~RESULTS/Methods/visual_17334_MW.pdf",
       plot = visual_17334_MW$plot,
       width = 40, height = 25, units = "cm")

# ------------------------------
#   computation time compariosn
# ------------------------------

time_17334_sample <- plot.time(data_17334, BIN=F, Shiftbin = F, rand = 1000)

time_17334_BIN <- plot.time(data_17334, BIN=T,BIN.size = 500, Shiftbin = F, rand = 1000)

time_17334_MW <- plot.time(data_17334, BIN=T,BIN.size = 500, Shiftbin = T, N.shifts = 5, rand = 1000)

time_17334_sum <- ggarrange(
  time_17334_sample$data %>%
    ggplot(aes(y=elapsed, x=smooth, fill=DC))+
    geom_hline(yintercept = c(50,100,150), color="gray80")+
    geom_bar(stat="identity", position = "dodge", color="gray30")+
    ggtitle(paste("N samples =",150,
                  ",N taxa =",80,
                  ",N.randomisation =",1000,
                  ",sample"))+
    theme_classic()+
    coord_cartesian(ylim = c(0,160))+
    ylab("computation time (s)"),
  time_17334_BIN$data%>%
    ggplot(aes(y=elapsed, x=smooth, fill=DC))+
    geom_hline(yintercept = c(50,100,150), color="gray80")+
    geom_bar(stat="identity", position = "dodge", color="gray30")+
    ggtitle(paste("N samples =",150,
                  ",N taxa =",80,
                  ",N.randomisation =",1000,
                  ",Binning"))+
    theme_classic()+
    coord_cartesian(ylim = c(0,160))+
    ylab("computation time (s)"),
  time_17334_MW$data%>%
    ggplot(aes(y=elapsed, x=smooth, fill=DC))+
    geom_hline(yintercept = c(50,100,150), color="gray80")+
    geom_bar(stat="identity", position = "dodge", color="gray30")+
    ggtitle(paste("N samples =",150,
                  ",N taxa =",80,
                  ",N.randomisation =",1000,
                  ",Moving window"))+
    theme_classic()+
    coord_cartesian(ylim = c(0,160))+
    ylab("computation time (s)"),
  nrow = 3, ncol = 1, common.legend = T, legend = "right"
)

time_17334_sum

ggsave("~RESULTS/Methods/time_17334_sum.pdf",
       plot = time_17334_sum,
       height = 18, width = 25, units="cm")



# ----------------------------------------------
#               SUMMARY 
# ----------------------------------------------

data_ld_sum <- rbind(
  data.frame(perform_sim_ld_recent_MW$data,Position="recent"),
  data.frame(perform_sim_ld_late_MW$data,Position="late")
) %>%
  filter(SIGNIF=="Peak.gam")

data_ld_sum <- within(data_ld_sum, SEGMENT <- factor(SEGMENT, levels = c("focus","empty")))

plot_sum_ld_MW_gam <- data_ld_sum %>%
  ggplot(aes(y=VALUE.M,x=Position,fill=SEGMENT, group=SEGMENT))+
  geom_bar(stat="identity", position="dodge", color="gray30")+
  geom_errorbar(aes(ymin=VALUE.05,ymax=VALUE.95, group=SEGMENT),
                position=position_dodge(width=0.9), width=0.2, size=0.5, color="gray30")+
  facet_grid(smooth~DC)+
  scale_fill_manual("Position in sequence", labels=c("focal area (correct detection)","outside of focal area (false positive)"),
                    values = c("darkseagreen","coral"))+
  ylab("Percentage of Peak detection")+xlab("Position of enviromental change")+
  theme_classic()+
  theme(legend.position = "right")
  
plot_sum_ld_MW_gam

ggsave("~RESULTS/Methods/plot_sum_ld_MW_gam.pdf",
       plot = plot_sum_ld_MW_gam,
       height = 25, width = 25, units="cm")


data_hd_sum <- rbind(
  data.frame(perform_sim_hd_recent_MW$data,Position="recent"),
  data.frame(perform_sim_hd_late_MW$data,Position="late")
) %>%
  filter(SIGNIF=="Peak.gam")

data_hd_sum <- within(data_hd_sum, SEGMENT <- factor(SEGMENT, levels = c("focus","empty")))

plot_sum_hd_MV_gam <- data_hd_sum %>%
  ggplot(aes(y=VALUE.M,x=Position,fill=SEGMENT, group=SEGMENT))+
  geom_bar(stat="identity", position="dodge", color="gray30")+
  geom_errorbar(aes(ymin=VALUE.05,ymax=VALUE.95, group=SEGMENT),
                position=position_dodge(width=0.9), width=0.2, size=0.5, color="gray30")+
  facet_grid(smooth~DC)+
  scale_fill_manual("Position in sequence", labels=c("focal area (correct detection)","outside of focal area (false positive)"),
                    values = c("darkseagreen","coral"))+
  ylab("Percentage of Peak detection")+xlab("Position of enviromental change")+
  theme_classic()+
  theme(legend.position = "right")


plot_sum_hd_MV_gam

ggsave("~RESULTS/Methods/plot_sum_hd_MV_gam.pdf",
       plot = plot_sum_hd_MV_gam,
       height = 25, width = 25, units="cm")

data_mag <- rbind(
  data.frame(visual_sim_ld_recent_MW$data, DIVERSITY="ld" ,POSITION = "recent"),
  data.frame(visual_sim_ld_late_MW$data, DIVERSITY="ld" ,POSITION = "late"),
  data.frame(visual_sim_hd_recent_MW$data, DIVERSITY="hd" ,POSITION = "recent"),
  data.frame(visual_sim_hd_late_MW$data, DIVERSITY="hd" ,POSITION = "late")
)

plot_mag_MW <- data_mag %>%
  select(RUN.RoC,smooth,DC,DIVERSITY,POSITION) %>%
  group_by(smooth,DC, DIVERSITY, POSITION) %>%
    summarise(RoC.max = max(RUN.RoC)) %>%
  ggplot(aes(y=RoC.max,x=POSITION,fill=DIVERSITY))+
  geom_hline(yintercept = c(0,1,2,10,20,30), color="gray80")+
  geom_bar(stat="identity", position = "dodge", color="gray30")+
  facet_grid(smooth~DC)+
  scale_y_continuous(trans = "log1p", breaks = c(0,1,2,10,20,30))+
  theme_classic()+
  scale_fill_manual("Diversity",labels=c("low","high"), values = c("#F8766D","#00BFC4"))+
  xlab("Position of enviromental change")+ylab("Maximum value of Rate of CHange")

plot_mag_MW

ggsave("~RESULTS/Methods/plot_mag_MW.pdf",
       plot = plot_mag_MW,
       height = 25, width = 25, units="cm")

# -----------------------------------------

# save.image("~/DATA/temp/ENV_METHOD_20200408.RData")
# load("~/DATA/temp/ENV_METHOD_20200408.RData")


# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())



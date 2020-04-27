# load("~/DATA/temp/ENV_METHOD_20200427.RData")

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


fc_calculate_RoC_comparison <- function(data, BIN, BIN.size, Shiftbin, N.shifts, rand, standardise ,interest.treshold)
{
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  
  for(i in 1:20)
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
    
    data.temp.sum<- data.frame(data.temp,performance.smooth[i],performance.DC[i])
    names(data.temp.sum) <- c(names(data.temp),"smooth","DC")
    
    if(i == 1){
      perfomance.tibble <- data.temp.sum} else {
        perfomance.tibble <- rbind(perfomance.tibble,data.temp.sum)
      }
  }
  return(perfomance.tibble)
}

fc_get_pollen_data <- function (data, sm.type, N.taxa)
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
    select(any_of(Common.list)) %>%
    rownames_to_column() %>%
    pivot_longer(cols = -c(rowname)) %>%
    rename(sample.id = rowname) %>%
    inner_join(.,data.ext$Age, by="sample.id")
  
  return (plot.data)
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
age_lim <- 8000


# -----------------------------------------
#
#             STATISTICAL COMP
# 
# -----------------------------------------

perform_sim_ld_recent_MW <- fc_test_simlutated_data_succsess(time= time_seq,
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


perform_sim_ld_late_MW <- fc_test_simlutated_data_succsess(time= time_seq,
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

perform_sim_hd_recent_MW <- fc_test_simlutated_data_succsess(time= time_seq,
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

perform_sim_hd_late_MW <- fc_test_simlutated_data_succsess(time= time_seq,
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



# -----------------------------------------
#
#       MAGNITUDE COMPARISON
# 
# -----------------------------------------

# low diversity - recent
mag_sim_ld_recent_MW <- fc_test_simlutated_data_magnitude(time= time_seq,
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
                                                          Shiftbin=T,
                                                          N.shifts=5,
                                                          rand.sets=N_rep,
                                                          interest.treshold=age_lim)

# low diversity - late
mag_sim_ld_late_MW <- fc_test_simlutated_data_magnitude(time= time_seq,
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

# high diversity - recent
mag_sim_hd_recent_MW <- fc_test_simlutated_data_magnitude(time= time_seq,
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
                                                          Shiftbin=T,
                                                          N.shifts=5,
                                                          rand.sets=N_rep,
                                                          interest.treshold=age_lim)

# high diversity - late
mag_sim_hd_late_MW <- fc_test_simlutated_data_magnitude(time= time_seq,
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




data_example <- list(dataset.id = tibble_Europe2$dataset.id[[2]],
                   filtered.counts = tibble_Europe2$filtered.counts[[2]],
                   list_ages = tibble_Europe2$list_ages[[2]])

data_example_MW <- fc_calculate_RoC_comparison(data_example,
                                   BIN = T, 
                                   BIN.size = 500,
                                   Shiftbin = T, 
                                   N.shifts = 5,
                                   rand = 1000,
                                   standardise = T,
                                   interest.treshold =  age_lim)




# ----------------------------------------------
#
#               FIGURES 
#
# ----------------------------------------------

##############
#   FIG 1
#############

FIG1_visual_example_MW <-data_example_MW %>%
  ggplot(aes(y=RUN.RoC, 
             x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0))+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = filter(data_example_MW, Peak.treshold.95==T ),color="blue", size=2, shape=1, alpha=2/3)+
  geom_point(data = filter(data_example_MW, Peak.gam==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_point(data = filter(data_example_MW, Peak.SNI==T ),color="red", size=2, shape=8, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  xlab("Age (cal yr BC)")+ylab("Rate of Change score")+
  facet_grid(smooth~DC, scales = "free_x")

FIG1_visual_example_MW

ggsave("~/RESULTS/Methods/FIN/FIG1_visual_example_MW.pdf",
       plot = FIG1_visual_example_MW,
       width = 20, height = 12, units = "cm")



##############
#   FIG 2
#############

data_mag_sum <- rbind(
  data.frame(mag_sim_ld_recent_MW$data,Position="recent",Diversity="low"),
  data.frame(mag_sim_ld_late_MW$data,Position="late",Diversity="low"),
  data.frame(mag_sim_hd_recent_MW$data,Position="recent",Diversity="high"),
  data.frame(mag_sim_hd_late_MW$data,Position="late",Diversity="high")
)

data_mag_sum$Dataset_type <- c(rep("LD-R",20),rep("LD-L",20),rep("HD-R",20),rep("HD-L",20))

data_mag_sum <- within(data_mag_sum, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
data_mag_sum <- within(data_mag_sum, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))
data_mag_sum <- within(data_mag_sum, Dataset_type <- factor(Dataset_type, levels = c("LD-R","LD-L","HD-R","HD-L")))

FIG2_mag_MW  <- data_mag_sum %>%
  ggplot(aes(y=RoC_max,fill=Dataset_type,x=Dataset_type))+
  geom_bar(stat="identity", position = "dodge", color="gray30")+
  geom_errorbar(aes(ymax=RoC_max+RoC_max_SD, ymin=RoC_max-RoC_max_SD), 
                position = position_dodge(width=0.9), width=0.2, size=0.5, color="gray50")+
  facet_grid(DC~SMOOTH, scales = "free_y")+
  theme_classic()+
  #scale_fill_manual("Position of enviromental change",labels=c("recent","late"), values = c("#F8766D","#00BFC4"))+
  #scale_x_discrete(label=c("low diversity","high diversity"))+
  theme(axis.text.x = element_text(angle = -45,hjust = 0.3, vjust = 0.2 ),
        legend.position = "none")+
  xlab("Type of simulated dataset")+ylab("Maximum value of Rate of Change")

FIG2_mag_MW

ggsave("~/RESULTS/Methods/FIN/FIG2_mag_MW.pdf",
       plot = FIG2_mag_MW,
       height = 12, width = 20, units="cm")

##############
#   FIG 3
#############


data_success_sum <- rbind(
  data.frame(perform_sim_ld_recent_MW$data,Position="recent", Diversity= "low"),
  data.frame(perform_sim_ld_late_MW$data,Position="late", Diversity= "low"),
  data.frame(perform_sim_hd_recent_MW$data,Position="recent", Diversity= "high"),
  data.frame(perform_sim_hd_late_MW$data,Position="late", Diversity= "high")
) %>%
  filter(SIGNIF=="Peak.gam")

data_success_sum$Dataset_type <- c(rep("LD-R",40),rep("LD-L",40),rep("HD-R",40),rep("HD-L",40))
data_success_sum <- within(data_success_sum, SEGMENT <- factor(SEGMENT, levels = c("focus","empty")))
data_success_sum <- within(data_success_sum, Dataset_type <- factor(Dataset_type, levels = c("LD-R","LD-L","HD-R","HD-L")))

FIG3_sum_MW_gam <- data_success_sum %>%
  ggplot(aes(y=VALUE.M,x=Dataset_type,fill=SEGMENT, group=SEGMENT))+
  geom_bar(stat="identity", position="dodge", color="gray30")+
  geom_errorbar(aes(ymin=VALUE.M-VALUE.SD,ymax=VALUE.M+VALUE.SD, group=SEGMENT),
                position=position_dodge(width=0.9), width=0.2, size=0.5, color="gray50")+
  facet_grid(smooth~DC)+
  scale_fill_manual("Position in sequence", labels=c("focal area (correct detection)","outside of focal area (false positive)"),
                    values = c("darkseagreen","coral"))+
  ylab("Percentage of Peak detection")+xlab("Type of simulated dataset")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45,hjust = 0.3, vjust = 0.2 ),
        legend.position = "bottom")

FIG3_sum_MW_gam

ggsave("~/RESULTS/Methods/FIN/FIG3_sum_MW_gam.pdf",
       plot = FIG3_sum_MW_gam,
       height = 12, width = 20, units="cm")

data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  filter(Position == "recent") %>%
  group_by(smooth,DC) %>%
  summarise(
    VALUE = mean(VALUE.M),
    SD = mean(VALUE.SD)
  ) %>%
  ungroup() %>%
  group_by(smooth) %>%
  arrange(-VALUE, .by_group=T) %>%
  View()


data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  filter(Position == "recent") %>%
  group_by(smooth) %>%
  summarise(
    VALUE = mean(VALUE.M),
    SD = mean(VALUE.SD)
  ) %>%
  View()


data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  filter(Position == "late") %>%
  group_by(smooth,DC) %>%
  summarise(
    VALUE = mean(VALUE.M)
  ) %>%
  ungroup() %>%
  group_by(smooth) %>%
  arrange(-VALUE, .by_group=T) %>%
  View()

data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  filter(Position == "late") %>%
  group_by(smooth,DC) %>%
  summarise(
    SD= mean(VALUE.SD)
  ) %>%
  ungroup() %>%
  group_by(smooth) %>%
  arrange(SD, .by_group=T) %>%
  View()


data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  filter(Position == "late") %>%
  group_by(smooth,DC) %>%
  summarise(
    VALUE = mean(VALUE.M),
    SD = mean(VALUE.SD)
  ) %>%
  ungroup() %>%
  group_by(smooth) %>%
  arrange(SD, .by_group=T) %>%
  View()

data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  filter(Position == "late") %>%
  summarise(
    VALUE = mean(VALUE.M),
    SD = mean(VALUE.SD)
  ) %>%
  View()

data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  filter(Position == "late") %>%
  group_by(smooth) %>%
  summarise(
    VALUE = mean(VALUE.M),
    SD = mean(VALUE.SD)
  ) %>%
  View()

data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  filter(Position == "late") %>%
  group_by(DC) %>%
  summarise(
    VALUE = mean(VALUE.M),
    SD = mean(VALUE.SD)
  ) %>%
  View()



data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  filter(DC == "chord") %>%
  group_by(Position, smooth) %>%
  summarise(
    VALUE = mean(VALUE.M),
    SD = sd(VALUE.M)
  ) %>%
  ungroup() %>%
  arrange(-VALUE) %>%
  View()

data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  filter(DC == "chord") %>%
  group_by(Position, smooth) %>%
  summarise(
    VALUE = mean(VALUE.M),
    SD = sd(VALUE.M)
  ) %>%
  ungroup() %>%
  arrange(VALUE) %>%
  View()




############
#   FIG 4
############


#SITE A
which(tibble_Europe2$dataset.id %in%  17334 )

data_site_A <- list(dataset.id = tibble_Europe2$dataset.id[[2]],
                    filtered.counts = tibble_Europe2$filtered.counts[[2]],
                    list_ages = tibble_Europe2$list_ages[[2]])

data_site_A_pollen <-fc_get_pollen_data(data_site_A, sm.type = "shep",N.taxa = 10)

data_site_A$filtered.counts %>%
  as_tibble() %>%
  filter(data_site_A$list_ages$ages$age < 9000) %>%
  dim()

data_site_A_pollen %>%
  group_by(name) %>%
  summarise(SUM = sum(value)) %>%
  arrange(-SUM)
  

#Site B

which(tibble_Europe2$dataset.id %in%  40951 )

data_site_B <- list(dataset.id = tibble_Europe2$dataset.id[[224]],
                    filtered.counts = tibble_Europe2$filtered.counts[[224]],
                    list_ages = tibble_Europe2$list_ages[[224]])


data_site_B_pollen <-fc_get_pollen_data(data_site_B, sm.type = "shep",N.taxa = 10)

data_site_B$filtered.counts %>%
  filter(data_site_B$list_ages$ages$age < 9000) %>%
  as_tibble() %>%
  dim()

data_site_B_pollen %>%
  group_by(name) %>%
  summarise(SUM = sum(value)) %>%
  arrange(-SUM)


# SITE c

which(tibble_Europe2$dataset.id %in%  4012 )

data_site_C <- list(dataset.id = tibble_Europe2$dataset.id[[45]],
                    filtered.counts = tibble_Europe2$filtered.counts[[45]],
                    list_ages = tibble_Europe2$list_ages[[45]])

data_site_C_pollen <-fc_get_pollen_data(data_site_C, sm.type = "shep",N.taxa = 10)

data_site_C$filtered.counts %>%
  as_tibble() %>%
  filter(data_site_C$list_ages$ages$age < 9000) %>%
  dim()

data_site_C_pollen %>%
  group_by(name) %>%
  summarise(SUM = sum(value)) %>%
  arrange(-SUM)


# Pollen taxa table

common_taxa<- c(data_site_A_pollen$name,
                data_site_B_pollen$name,
                data_site_C_pollen$name) %>%
  unique()

library (RColorBrewer)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
Palette.1<- getPalette(length(common_taxa))
names(Palette.1)<- sort(common_taxa)  


# Rate of Change

data_site_A_RoC <- fc_ratepol(data.source.pollen = data_site_A$filtered.counts,
                              data.source.age = data_site_A$list_ages,
                              sm.type = "shep", 
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              BIN = T,
                              BIN.size = 500,
                              Shiftbin  = T,
                              N.shifts = 5,
                              rand = 1000,
                              standardise = T, 
                              S.value = 150, 
                              DC = "chord",
                              interest.treshold = age_lim,
                              Debug = F)


data_Site_C_RoC <- fc_ratepol(data.source.pollen = data_site_C$filtered.counts,
                              data.source.age = data_site_C$list_ages,
                              sm.type = "shep", 
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              BIN = T,
                              BIN.size = 500,
                              Shiftbin  = T,
                              N.shifts = 5,
                              rand = 1000,
                              standardise = T, 
                              S.value = 150, 
                              DC = "chord",
                              interest.treshold = age_lim,
                              Debug = F)


data_Site_B_RoC <- fc_ratepol(data.source.pollen = data_site_B$filtered.counts,
                              data.source.age = data_site_B$list_ages,
                              sm.type = "shep", 
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              BIN = T,
                              BIN.size = 500,
                              Shiftbin  = T,
                              N.shifts = 5,
                              rand = 1000,
                              standardise = T, 
                              S.value = 150, 
                              DC = "chord",
                              interest.treshold = age_lim,
                              Debug = F)


# FIGURES 

FIG4_Site_A_1 <- data_site_A$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = seq(from=0,to=3e-4, by=1e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = seq(from=0,to=3e-4, by=1e-4))+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))+
  xlab("Age (cal yr BC)")+ylab("")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


FIG4_Site_B_1 <- data_site_B$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = seq(from=0,to=3e-4, by=1e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = seq(from=0,to=3e-4, by=1e-4))+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))+
  xlab("Age (cal yr BC)")+ylab("")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


FIG4_Site_C_1 <- data_site_C$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = seq(from=0,to=3e-4, by=1e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))+
  scale_y_continuous(breaks = seq(from=0,to=3e-4, by=1e-4))+
  xlab("Age (cal yr BC)")+ylab("Density of samples")


library(cowplot)
my_legend <- cowplot::get_legend(ggplot(data = data.frame(NAME=common_taxa,X=1), aes(x=X, fill=NAME))+
                             geom_density(alpha=1/3)+
                             theme_classic()+
                             theme(legend.position = "bottom")+
                             scale_fill_manual("pollen taxa",values = Palette.1))

plot(my_legend)

FIG4_Site_A_2 <-  data_site_A_pollen %>%
  ggplot(aes( y=value, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_hline(yintercept = seq(from=0,to=1, by=0.25), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
              color="gray20", alpha=1/5, size=0.1)+
  scale_fill_manual("pollen taxa",values = Palette.1)+
  xlab("")+ylab("")+
  coord_flip(xlim=c(age_lim,0), ylim = c(0,1))+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

FIG4_Site_B_2 <-  data_site_B_pollen %>%
  ggplot(aes( y=value, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_hline(yintercept = seq(from=0,to=1, by=0.25), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
              color="gray20", alpha=1/5, size=0.1)+
  scale_fill_manual("pollen taxa",values = Palette.1, drop=FALSE)+
  xlab("")+ylab("")+
  coord_flip(xlim=c(age_lim,0), ylim = c(0,1))+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

FIG4_Site_C_2 <-  data_site_C_pollen %>%
  ggplot(aes( y=value, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_hline(yintercept = seq(from=0,to=1, by=0.25), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
              color="gray20", alpha=1/5, size=0.1)+
  scale_fill_manual("pollen taxa",values = Palette.1)+
  xlab("")+ylab("% of pollen grains")+
  coord_flip(xlim=c(age_lim,0), ylim = c(0,1))+
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


FIG4_Site_A_3 <- data_site_A_RoC %>%
ggplot(aes(y=RUN.RoC, 
           x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,1.5))+
  geom_hline(yintercept = seq(from=0,to=1.5, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = filter(data_site_A_RoC, Peak.gam==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  xlab("")+ylab("")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank())


FIG4_Site_B_3 <- data_Site_B_RoC %>%
  ggplot(aes(y=RUN.RoC, 
             x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,1.5))+
  geom_hline(yintercept = seq(from=0,to=1.5, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = filter(data_Site_B_RoC, Peak.gam==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  xlab("")+ylab("")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank())


FIG4_Site_C_3 <- data_Site_C_RoC %>%
  ggplot(aes(y=RUN.RoC, 
             x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,1.5))+
  geom_hline(yintercept = seq(from=0,to=1.5, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = filter(data_Site_C_RoC, Peak.gam==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  xlab("")+ylab("Rate of Change score")+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


FIG4_Site_comparison <- ggarrange(
  FIG4_Site_A_1,FIG4_Site_A_2,FIG4_Site_A_3,
  FIG4_Site_B_1,FIG4_Site_B_2,FIG4_Site_B_3,
  FIG4_Site_C_1,FIG4_Site_C_2,FIG4_Site_C_3,
  #labels = c(17334,"","",40951,"","",4012,"",""),
  labels = c("A","","","B","","","C","",""),
  ncol = 3, nrow = 3, legend = "none", align = "hv")


FIG4_Site_comparison_legend <- ggarrange(
  FIG4_Site_comparison,
  my_legend, 
  ncol=1, heights = c(10,1))

ggsave("~/RESULTS/Methods/FIN/FIG4_Site_comparison.pdf",
       plot = FIG4_Site_comparison_legend, 
       height = 20, width = 20, units="cm")



###############
# SUPLEMENTARY 
###############

# Need to go to the DEBUG and run SIMULATION MANUALY!!!


Supplementary_F1 <-ggarrange(as.data.frame(forcing) %>%
                               mutate(AGE = time) %>%
                               pivot_longer(., cols = -c(AGE)) %>%
                               arrange(AGE,value) %>%
                               ggplot(aes(x=AGE, y= value))+
                               geom_vline(xintercept = breaks, color="gray80", size=0.1)+
                               geom_line(aes(color=name))+
                               theme_classic()+
                               coord_flip(xlim=c(8000,0))+
                               scale_x_continuous(trans = "reverse")+
                               theme(legend.position = "none")+
                               xlab("Age (cal yr BC)"),
                             fc_extract(data.source.pollen, data.source.age) %>%
                               fc_smooth("none") %>%
                               fc_check(., proportion = T) %>%
                               pluck("Pollen") %>%
                               mutate(AGE = time) %>%
                               pivot_longer(-c(AGE)) %>%
                               ggplot(aes(x=AGE, y=value))+
                               geom_ribbon(aes(ymin=rep(0,length(value)), ymax=value, fill=name), 
                                           color="gray20", alpha=1/5, size=0.1)+
                               geom_vline(xintercept = breaks, color="gray80", size=0.1)+
                               coord_flip(xlim=c(8000,0), ylim=c(0,1))+
                               theme_classic()+
                               scale_x_continuous(trans = "reverse")+
                               ylab("% of pollen grains")+xlab("")+
                               theme(legend.position = "none",
                                     axis.ticks.y = element_blank(),
                                     axis.text.y = element_blank()),
                             ncol=2, align = "h")


Supplementary_F1

ggsave("~/RESULTS/Methods/FIN/Supplementary_F1.pdf",
       plot = Supplementary_F1,
       height = 10, width = 15, units="cm")




data_supp_fig_MW <- rbind(
  perform_sim_ld_recent_MW$data,
  perform_sim_ld_late_MW$data,
  perform_sim_hd_recent_MW$data,
  perform_sim_hd_late_MW$data
)

data_supp_fig_MW$Dataset_type <- c(rep("LD-R",120),rep("LD-L",120),rep("HD-R",120),rep("HD-L",120))

data_supp_fig_MW <- within(data_supp_fig_MW, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
data_supp_fig_MW <- within(data_supp_fig_MW, smooth <- factor(smooth, levels = c("none","m.avg","grim","age.w","shep")))
data_supp_fig_MW <- within(data_supp_fig_MW, Dataset_type <- factor(Dataset_type, levels = c("LD-R","LD-L","HD-R","HD-L")))
data_supp_fig_MW <- within(data_supp_fig_MW, SEGMENT <- factor(SEGMENT, levels = c("focus","empty")))
data_supp_fig_MW <- within(data_supp_fig_MW, SIGNIF <- factor(SIGNIF, levels = c("Peak.treshold","Peak.gam","Peak.SNI")))
levels(data_supp_fig_MW$SIGNIF) <- c("Threshold","GAM","SNI")


Supplementary_F2 <- data_supp_fig_MW %>%
  ggplot(aes(y=VALUE.M,x=Dataset_type,fill=SEGMENT, group=SEGMENT))+
  geom_bar(stat="identity", position="dodge", color="gray30")+
  geom_errorbar(aes(ymin=VALUE.M-VALUE.SD,ymax=VALUE.M+VALUE.SD, group=SEGMENT),
                position=position_dodge(width=0.9), width=0.2, size=0.5, color="gray50")+
  facet_grid(smooth~DC+SIGNIF)+
  scale_fill_manual("Position in sequence", labels=c("focal area (correct detection)","outside of focal area (false positive)"),
                    values = c("darkseagreen","coral"))+
  ylab("Percentage of Peak detection")+xlab("Type of simulated dataset")+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = -45,hjust = 0.3, vjust = 0.2 ))


Supplementary_F2

ggsave("~/RESULTS/Methods/FIN/Supplementary_F2.pdf",
       plot = Supplementary_F1,
       height = 15, width = 22, units="cm")


# save.image("~/DATA/temp/ENV_METHOD_20200427.RData")

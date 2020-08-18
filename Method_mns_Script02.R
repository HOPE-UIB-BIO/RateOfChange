# load("C:/Users/omo084/OneDrive - University of Bergen/PRIVATE/ROC_method/ENV_METHOD_2020807.RData")

# ----------------------------------------------
#                     SETUP
# ----------------------------------------------
library(tidyverse)
library(reshape2)
library(devtools)
library(ggpubr)
library(doSNOW)
library(parallel)
library(foreach)
library(doParallel)
library(scales)
library(mgcv)
library(maps)
library(RColorBrewer)
library(MuMIn)
library(glmmTMB)
library(emmeans)

theme_set(theme_classic())

# ----------------------------------------------
#                 LOAD DATA 
# ----------------------------------------------

Hope_smooth <- readRDS("~/DATA/input/HOPE_data_smooth20200728.RDS")

files.sources <- list.files("~/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

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

# Number of simulated datasest of pollen data
N_rep <- 100

# template of time sequence with uneven distribution of points
time_seq <- Hope_smooth$list_ages[[413]]$ages$age

# maximal time of focus 
age_lim <- 8000


# -----------------------------------------
#
#             STATISTICAL COMP
# 
# -----------------------------------------

# data simulation

sim_ld_recent <- fc_simulate_pollen_data_in_multiple_datasets(time=time_seq, 
                                                              nforc=N_env, 
                                                              mean=100, 
                                                              sdev=.15, 
                                                              nprox=low_diversity, 
                                                              var=20, 
                                                              range=15,
                                                              manual.edit = T,
                                                              breaks=breaks_recent,
                                                              jitter = T,
                                                              rarity=T,
                                                              N.datasets=N_rep)

sim_ld_recent_levels <- fc_simulate_pollen_data_in_all_methods(sim_ld_recent,
                                                           Working.Unit="levels", 
                                                           interest.treshold=8000)

sim_ld_recent_BINs <- fc_simulate_pollen_data_in_all_methods(sim_ld_recent,
                                                           Working.Unit="BINs", 
                                                           BIN.size=500,
                                                           interest.treshold=8000)

sim_ld_recent_MW <- fc_simulate_pollen_data_in_all_methods(sim_ld_recent,
                                                           Working.Unit="MW", 
                                                           BIN.size=500, 
                                                           N.shifts=5,
                                                           interest.treshold=8000)


sim_ld_late <- fc_simulate_pollen_data_in_multiple_datasets(time=time_seq, 
                                                              nforc=N_env, 
                                                              mean=100, 
                                                              sdev=.15, 
                                                              nprox=high_diversity, 
                                                              var=20, 
                                                              range=15,
                                                              manual.edit = T,
                                                              breaks=breaks_late,
                                                              jitter = T,
                                                              rarity=T,
                                                              N.datasets=N_rep)

sim_ld_late_levels <- fc_simulate_pollen_data_in_all_methods(sim_ld_late,
                                                               Working.Unit="levels", 
                                                               interest.treshold=8000)

sim_ld_late_BINs <- fc_simulate_pollen_data_in_all_methods(sim_ld_late,
                                                             Working.Unit="BINs", 
                                                             BIN.size=500,
                                                             interest.treshold=8000)

sim_ld_late_MW <- fc_simulate_pollen_data_in_all_methods(sim_ld_late,
                                                           Working.Unit="MW", 
                                                           BIN.size=500, 
                                                           N.shifts=5,
                                                           interest.treshold=8000)

sim_hd_recent <- fc_simulate_pollen_data_in_multiple_datasets(time=time_seq, 
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
                                                              N.datasets=N_rep)

sim_hd_recent_levels <- fc_simulate_pollen_data_in_all_methods(sim_hd_recent,
                                                               Working.Unit="levels", 
                                                               interest.treshold=8000)

sim_hd_recent_BINs <- fc_simulate_pollen_data_in_all_methods(sim_hd_recent,
                                                             Working.Unit="BINs", 
                                                             BIN.size=500,
                                                             interest.treshold=8000)

sim_hd_recent_MW <- fc_simulate_pollen_data_in_all_methods(sim_hd_recent,
                                                           Working.Unit="MW", 
                                                           BIN.size=500, 
                                                           N.shifts=5,
                                                           interest.treshold=8000)


sim_hd_late <- fc_simulate_pollen_data_in_multiple_datasets(time=time_seq, 
                                                            nforc=N_env, 
                                                            mean=100, 
                                                            sdev=.15, 
                                                            nprox=high_diversity, 
                                                            var=20, 
                                                            range=15,
                                                            manual.edit = T,
                                                            breaks=breaks_late,
                                                            jitter = T,
                                                            rarity=T,
                                                            N.datasets=N_rep)

sim_hd_late_levels <- fc_simulate_pollen_data_in_all_methods(sim_hd_late,
                                                             Working.Unit="levels", 
                                                             interest.treshold=8000)

sim_hd_late_BINs <- fc_simulate_pollen_data_in_all_methods(sim_hd_late,
                                                           Working.Unit="BINs", 
                                                           BIN.size=500,
                                                           interest.treshold=8000)

sim_hd_late_MW <- fc_simulate_pollen_data_in_all_methods(sim_hd_late,
                                                         Working.Unit="MW", 
                                                         BIN.size=500, 
                                                         N.shifts=5,
                                                         interest.treshold=8000)

# -----------------------------------------
#
#       SUCCESS COMPARISON
# 
# -----------------------------------------

perform_sim_ld_recent_MW <- fc_test_simlutated_data_succsess(sim_ld_recent_MW, breaks = breaks_recent)
perform_sim_ld_late_MW <- fc_test_simlutated_data_succsess(sim_ld_late_MW, breaks = breaks_late)
perform_sim_hd_recent_MW <- fc_test_simlutated_data_succsess(sim_hd_recent_MW, breaks = breaks_recent)
perform_sim_hd_late_MW <- fc_test_simlutated_data_succsess(sim_hd_late_MW, breaks = breaks_late)

perform_sim_ld_recent_BINs <- fc_test_simlutated_data_succsess(sim_ld_recent_BINs, breaks = breaks_recent)
perform_sim_ld_late_BINs <- fc_test_simlutated_data_succsess(sim_ld_late_BINs, breaks = breaks_late)
perform_sim_hd_recent_BINs <- fc_test_simlutated_data_succsess(sim_hd_recent_BINs, breaks = breaks_recent)
perform_sim_hd_late_BINs <- fc_test_simlutated_data_succsess(sim_hd_late_BINs, breaks = breaks_late)

perform_sim_ld_recent_levels <- fc_test_simlutated_data_succsess(sim_ld_recent_levels, breaks = breaks_recent)
perform_sim_ld_late_levels <- fc_test_simlutated_data_succsess(sim_ld_late_levels, breaks = breaks_late)
perform_sim_hd_recent_levels <- fc_test_simlutated_data_succsess(sim_hd_recent_levels, breaks = breaks_recent)
perform_sim_hd_late_levels <- fc_test_simlutated_data_succsess(sim_hd_late_levels, breaks = breaks_late)

data_success_sum <- rbind(
  data.frame(perform_sim_ld_recent_MW$RawData,Position="recent", Diversity= "low", WU="MW"),
  data.frame(perform_sim_ld_late_MW$RawData,Position="late", Diversity= "low", WU="MW"),
  data.frame(perform_sim_hd_recent_MW$RawData,Position="recent", Diversity= "high", WU="MW"),
  data.frame(perform_sim_hd_late_MW$RawData,Position="late", Diversity= "high", WU="MW"),
  data.frame(perform_sim_ld_recent_BINs$RawData,Position="recent", Diversity= "low", WU="BINs"),
  data.frame(perform_sim_ld_late_BINs$RawData,Position="late", Diversity= "low", WU="BINs"),
  data.frame(perform_sim_hd_recent_BINs$RawData,Position="recent", Diversity= "high", WU="BINs"),
  data.frame(perform_sim_hd_late_BINs$RawData,Position="late", Diversity= "high", WU="BINs"),
  data.frame(perform_sim_ld_recent_levels$RawData,Position="recent", Diversity= "low", WU="levels"),
  data.frame(perform_sim_ld_late_levels$RawData,Position="late", Diversity= "low", WU="levels"),
  data.frame(perform_sim_hd_recent_levels$RawData,Position="recent", Diversity= "high", WU="levels"),
  data.frame(perform_sim_hd_late_levels$RawData,Position="late", Diversity= "high", WU="levels")
) %>%
  as_tibble()

data_success_sum <- within(data_success_sum, Position <- factor(Position, levels = c("recent","late")))
levels(data_success_sum$Position) <- c("high density level","low density level")

data_success_sum <- within(data_success_sum, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))
levels(data_success_sum$SMOOTH) <- c("None","M.avg","Grimm","Age.w","Shep")

data_success_sum <- within(data_success_sum, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
levels(data_success_sum$DC) <- c("Euc","Euc.sd","Chord","Chisq")

data_success_sum <- within(data_success_sum, Diversity <- factor(Diversity, levels = c("low","high")))
levels(data_success_sum$Diversity) <- c("low richness","high richness")

data_success_sum <- within(data_success_sum, WU <- factor(WU, levels = c("levels","BINs","MW")))


# Produce FIG S2 !!!

# cluster setup

nrCores <- detectCores()

## FOCUS


data_success_focus <- data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  ungroup() %>%
  dplyr::select(-c(SEGMENT))

mod_success_focus <-  glmmTMB(VALUE~WU*PEAK*Position*Diversity+(WU|dataset.ID),
                           data=data_success_focus,
                           family=betabinomial(link = "logit"))




cl <- makeCluster(nrCores-1)
registerDoParallel(cl); 
clusterExport(cl,c("data_success_focus","nrCores"),envir=environment());
clusterEvalQ(cl,library("glmmTMB"))

mod_success_focus_dd <- pdredge(mod_success_focus,cluster = cl, trace = T)

(mod_success_focus_dd$AICc - mod_success_focus_dd$AICc[1])<2

mod_success_focus_dd[1,]
#  -> FULL model 

# Drop concirm the same thing
drop1(mod_success_focus, test = "Chisq", trace = T)


mod_success_focus_dd %>%
  as_tibble() %>%
  write.csv(.,"METHOD_RESULTS/mod_success_focus.csv")


data_success_focus_MW_G <- data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  filter(WU == "MW") %>% 
  filter(PEAK == "PEAK.G") %>%
  ungroup() %>%
  dplyr::select(-c(SEGMENT,WU,PEAK))

mod_success_focus_MW_G <-  glmmTMB(VALUE~Position*Diversity*DC*SMOOTH+(1|dataset.ID),
                                   data=data_success_focus_MW_G,
                                   family=betabinomial(link = "logit"))

cl <- makeCluster(nrCores-1)
registerDoParallel(cl); 
clusterExport(cl,c("data_success_focus_MW_G","nrCores"),envir=environment());
clusterEvalQ(cl,library("glmmTMB"))

mod_success_focus_MW_G_dd <- pdredge(mod_success_focus_MW_G,cluster = cl, trace = T)

(mod_success_focus_MW_G_dd$AICc - mod_success_focus_MW_G_dd$AICc[1])<2 

importance(mod_success_focus_MW_G_dd[1:5,])

mod_success_focus_MW_G_dd %>%
  as_tibble() %>%
  write.csv(.,"METHOD_RESULTS/mod_success_focus_MW_G_dd.csv")
  

# subset
mod_success_focus_MW_G_sub <- glmmTMB(VALUE~Position+Diversity+DC+SMOOTH+
                                        Diversity:Position+Diversity:SMOOTH+Position:SMOOTH+Diversity:Position:SMOOTH+DC:Position++DC:Diversity+
                                        (1|dataset.ID),
                                      data=data_success_focus_MW_G,
                                      family=betabinomial(link = "logit"),
                                      control = glmmTMBControl(parallel = nrCores-1,
                                                               optArgs = list(method="BFGS"),
                                                               optimizer = optim))


## EMPTY 

data_success_false <- data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  ungroup() %>%
  dplyr::select(-c(SEGMENT))

mod_success_false <-  glmmTMB((1-VALUE)~WU*PEAK*Position*Diversity+(WU|dataset.ID),
                              data=data_success_false,
                              family=betabinomial(link = "logit"))

cl <- makeCluster(nrCores-1)
registerDoParallel(cl); 
clusterExport(cl,c("data_success_false","nrCores"),envir=environment());
clusterEvalQ(cl,library("glmmTMB"))

mod_success_false_dd <- pdredge(mod_success_false,cluster = cl, trace = T)

(mod_success_false_dd$AICc - mod_success_false_dd$AICc[1])<2

importance(mod_success_false_dd[1:3,])


mod_success_false_dd %>%
  as_tibble() %>%
  write.csv(.,"METHOD_RESULTS/mod_success_false_dd.csv")


# sub
mod_success_false_sub <-  glmmTMB(VALUE~PEAK+Position+WU+Diversity+
                                    PEAK:Position+PEAK:WU+Position:WU+PEAK:Position:WU+
                                    (1|dataset.ID),
                              data=data_success_false,
                              family=betabinomial(link = "logit"),
                              control = glmmTMBControl(parallel = nrCores-1))

data_success_false_MW_G <- data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  filter(WU == "MW") %>% 
  filter(PEAK == "PEAK.G") %>%
  ungroup() %>%
  dplyr::select(-c(SEGMENT,WU,PEAK))

mod_success_false_MW_G <-  glmmTMB(VALUE~Position*Diversity*DC*SMOOTH+(1|dataset.ID),
                                   data=data_success_false_MW_G,
                                   family=betabinomial(link = "logit"))


cl <- makeCluster(nrCores-1)
registerDoParallel(cl); 
clusterExport(cl,c("data_success_false_MW_G","nrCores"),envir=environment());
clusterEvalQ(cl,library("glmmTMB"))

mod_success_false_MW_G_dd <- pdredge(mod_success_false_MW_G,cluster = cl, trace = T)

(mod_success_false_MW_G_dd$AICc - mod_success_false_MW_G_dd$AICc[1])<4 

importance(mod_success_false_MW_G_dd[1:6,])


mod_success_false_MW_G_dd %>%
  as_tibble() %>%
  write.csv(.,"METHOD_RESULTS/mod_success_false_MW_G_dd.csv")


# sub
mod_success_false_MW_G_sub <- glmmTMB(VALUE~Position+SMOOTH+Diversity+DC+Position:SMOOTH+Diversity:Position+
                                        (1|dataset.ID),
                                      data=data_success_false_MW_G,
                                      family=betabinomial(link = "logit"),
                                      control = glmmTMBControl(parallel = nrCores-1))

#optArgs = list(method="BFGS"),optimizer = optim


# -----------------------------------------
#
#       MAGNITUDE COMPARISON
# 
# -----------------------------------------

mag_sim_ld_recent_MW <- fc_test_simlutated_data_magnitude(sim_ld_recent_MW)
mag_sim_ld_late_MW <- fc_test_simlutated_data_magnitude(sim_ld_late_MW)
mag_sim_hd_recent_MW <- fc_test_simlutated_data_magnitude(sim_hd_recent_MW)
mag_sim_hd_late_MW <- fc_test_simlutated_data_magnitude(sim_hd_late_MW)

mag_sim_ld_recent_BINs <- fc_test_simlutated_data_magnitude(sim_ld_recent_BINs)
mag_sim_ld_late_BINs <- fc_test_simlutated_data_magnitude(sim_ld_late_BINs)
mag_sim_hd_recent_BINs <- fc_test_simlutated_data_magnitude(sim_hd_recent_BINs)
mag_sim_hd_late_BINs <- fc_test_simlutated_data_magnitude(sim_hd_late_BINs)

mag_sim_ld_recent_levels <- fc_test_simlutated_data_magnitude(sim_ld_recent_levels)
mag_sim_ld_late_levels <- fc_test_simlutated_data_magnitude(sim_ld_late_levels)
mag_sim_hd_recent_levels <- fc_test_simlutated_data_magnitude(sim_hd_recent_levels)
mag_sim_hd_late_levels <- fc_test_simlutated_data_magnitude(sim_hd_late_levels)


data_mag_sum <- rbind(
  data.frame(mag_sim_ld_recent_MW,Position="recent",Diversity="low",WU = "MW"),
  data.frame(mag_sim_ld_late_MW,Position="late",Diversity="low",WU = "MW"),
  data.frame(mag_sim_hd_recent_MW,Position="recent",Diversity="high",WU = "MW"),
  data.frame(mag_sim_hd_late_MW,Position="late",Diversity="high",WU = "MW"),
  data.frame(mag_sim_ld_recent_BINs,Position="recent",Diversity="low",WU = "BINs"),
  data.frame(mag_sim_ld_late_BINs,Position="late",Diversity="low",WU = "BINs"),
  data.frame(mag_sim_hd_recent_BINs,Position="recent",Diversity="high",WU = "BINs"),
  data.frame(mag_sim_hd_late_BINs,Position="late",Diversity="high",WU = "BINs"),
  data.frame(mag_sim_ld_recent_levels,Position="recent",Diversity="low",WU = "levels"),
  data.frame(mag_sim_ld_late_levels,Position="late",Diversity="low",WU = "levels"),
  data.frame(mag_sim_hd_recent_levels,Position="recent",Diversity="high",WU = "levels"),
  data.frame(mag_sim_hd_late_levels,Position="late",Diversity="high",WU = "levels")
) %>% as_tibble()

data_mag_sum <- within(data_mag_sum, Position <- factor(Position, levels = c("recent","late")))
levels(data_mag_sum$Position) <- c("high density level","low density level")

data_mag_sum <- within(data_mag_sum, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))
levels(data_mag_sum$SMOOTH) <- c("None","M.avg","Grimm","Age.w","Shep")

data_mag_sum <- within(data_mag_sum, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
levels(data_mag_sum$DC) <- c("Euc","Euc.sd","Chord","Chisq")

data_mag_sum <- within(data_mag_sum, Diversity <- factor(Diversity, levels = c("low","high")))
levels(data_mag_sum$Diversity) <- c("low richness","high richness")

cl <- makeCluster(nrCores-1)
registerDoParallel(cl); 
clusterExport(cl,c("data_mag_sum","nrCores"),envir=environment());
clusterEvalQ(cl,library("glmmTMB"))

# upper quantile
mod_mag_MW_upq <-  glmmTMB(RoC_upq~WU+Position+Diversity+DC+SMOOTH+ #5
                             WU:Position+WU:Diversity+WU:DC+WU:SMOOTH+ #4
                             Position:Diversity+Position:DC+Position:SMOOTH+ #3
                             Diversity:DC+Diversity:SMOOTH+ #2
                             DC:SMOOTH+ #1
                             WU:Position:Diversity+WU:Position:DC+WU:Position:SMOOTH+WU:Diversity:DC+WU:Diversity:SMOOTH+WU:DC:SMOOTH+
                             Position:Diversity:DC+Position:Diversity:SMOOTH+Position:DC:SMOOTH+Diversity:DC:SMOOTH+
                             (WU|dataset.ID),
                           data=data_mag_sum,
                           family=Gamma(link = "inverse"))

mod_mag_MW_upq_dd <- pdredge(mod_mag_MW_upq, trace = T, cluster = cl)

# median
mod_mag_MW_median <-  glmmTMB(RoC_median~WU+Position+Diversity+DC+SMOOTH+ #5
                                WU:Position+WU:Diversity+WU:DC+WU:SMOOTH+ #4
                                Position:Diversity+Position:DC+Position:SMOOTH+ #3
                                Diversity:DC+Diversity:SMOOTH+ #2
                                DC:SMOOTH+ #1
                                WU:Position:Diversity+WU:Position:DC+WU:Position:SMOOTH+WU:Diversity:DC+WU:Diversity:SMOOTH+WU:DC:SMOOTH+
                                Position:Diversity:DC+Position:Diversity:SMOOTH+Position:DC:SMOOTH+Diversity:DC:SMOOTH+
                                (WU|dataset.ID),
                              data=data_mag_sum,
                              family=Gamma(link = "inverse"))

mod_mag_MW_median_dd <- pdredge(mod_mag_MW_median, trace = T, cluster = cl)


data_mag_sum_MW <- data_mag_sum %>%
  filter(WU =="MW")

# max
mod_mag_MW_max <-  glmmTMB(RoC_max~Position*Diversity*DC*SMOOTH+ (1|dataset.ID),
                           data=data_mag_sum_MW,
                           family=Gamma(link = "inverse"))

mod_mag_MW_max_dd <- pdredge(mod_mag_MW_max, trace = T, cluster = cl)


emmeans(mod_mag_MW_upq, ~ WU+Position+Diversity+DC+SMOOTH,
        type = "response") %>%
  as_tibble() %>%
  mutate(PD = paste(Position,Diversity, sep = " - ")) %>%
  mutate(WU = factor(WU, levels = c("levels","BINs","MW"))) %>%
  mutate(WU = fct_recode(WU, "Mowing window" = "MW")) %>% 
  ggplot(aes(y=response,x=PD, color=SMOOTH,fill=SMOOTH, ymin=lower.CL, ymax=upper.CL))+
  facet_grid(WU~DC)+
  geom_hline(yintercept = seq(0,1.5,0.5), color="gray90")+
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5)+
  geom_errorbar(width=0.2, position = position_dodge(width = 0.5))+
  geom_point(shape=15,position = position_dodge(width = 0.5))+
  #coord_cartesian(ylim = c(0,1))+
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05, vjust = 1))+
  labs(y="Rate of change score",
       x= "",
       fill="",
       color="")


emmeans(mod_mag_MW_max, ~ Position+Diversity+DC+SMOOTH,
        type = "response") %>%
  as_tibble() %>%
  mutate(PD = paste(Position,Diversity, sep = " - ")) %>%
  ggplot(aes(y=response,x=DC, color=SMOOTH,fill=SMOOTH, ymin=lower.CL, ymax=upper.CL))+
  facet_grid(~PD)+
  geom_hline(yintercept = seq(0,2,0.5), color="gray90")+
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5)+
  geom_errorbar(width=0.2, position = position_dodge(width = 0.5))+
  geom_point(shape=15,position = position_dodge(width = 0.5))+
  #coord_cartesian(ylim = c(0,1))+
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05, vjust = 1))+
  labs(y="Rate of change score",
       x= "",
       fill="",
       color="")
# ----------------------------------------------
#
#               FIGURES 
#
# ----------------------------------------------


Color.legen_smooth <- brewer.pal(n = 5, name = 'Set2')
names(Color.legen_smooth) <- c("None","M.avg","Grimm","Age.w","Shep")


Color.legen_DC <- brewer.pal(n =4, name = 'Set3')
names(Color.legen_DC) <- c("Euc","Euc.sd","Chord","Chisq")

Color.legen_Position <- brewer.pal(n = 2, name = 'Set1')
names(Color.legen_Position) <- c("high density level","low density level")

Color.legen_Diversity <- brewer.pal(n = 2, name = 'Paired')
names(Color.legen_Diversity) <- c("low richness","high richness")





##############
#   FIG 1
#############


##############
#   FIG 2
#############

ggarrange(rbind(tibble(emmeans(mod_success_focus, ~ WU*PEAK*Position*Diversity,
                               type = "response") %>%
                         as_tibble(),METHOD = "CORRECT"),
                tibble(emmeans(mod_success_false_sub, ~ WU+PEAK+Position+Diversity+
                                 PEAK:Position+PEAK:WU+Position:WU+PEAK:Position:WU+Diversity,
                               type = "response") %>%
                         as_tibble(), METHOD = "FALSE")
) %>%
  mutate(PEAK = factor(PEAK, levels = c("PEAK.T","PEAK.G","PEAK.S"))) %>%
  mutate(PEAK = fct_recode(PEAK, 
                           "Threshold" = "PEAK.T",
                           "GAM" = "PEAK.G",
                           "SNI" = "PEAK.S")) %>%
  mutate(PD = paste(Position,Diversity, sep = " - ")) %>%
  mutate(PD = as.factor(PD)) %>%
  mutate(PD = fct_recode(PD,
                         "HR - R" = "high density level - high richness",
                         "LR - R" = "high density level - low richness",
                         "HR - L" = "low density level - high richness",
                         "LR - L" = "low density level - low richness",)) %>%
  ggplot(aes(y=prob,x=PD,color=METHOD, fill=METHOD, ymin=lower.CL, ymax=upper.CL))+
  facet_grid(WU~PEAK)+
  geom_hline(yintercept = seq(0,1,0.25), color="gray90",size=0.1)+
  geom_errorbar(width=0.2,color="gray60", position = position_dodge(width = 0.5),size=0.1)+
  geom_bar(stat = "identity",color="gray60", orientation="x", width = 0.5, position = position_dodge(width = 0.5),size=0.1)+
  #geom_point(shape=15, position = position_dodge(width = 0.5))+
  scale_fill_manual("Position in sequence", labels=c("Focal area (correct detection)","Outside of focal area (false positive)"),
                    values = c("darkseagreen","coral"))+
  coord_cartesian(ylim = c(0,1))+
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05, vjust = 1))+
  labs(y="",
       x= ""),
rbind(tibble(emmeans(mod_success_focus_MW_G_sub, ~ Position+Diversity+DC+SMOOTH+
                       Diversity:Position+Diversity:SMOOTH+Position:SMOOTH+Diversity:Position:SMOOTH+DC:Position++DC:Diversity,
                     type = "response") %>%
               as_tibble(), METHOD = "CORRECT"),
      tibble(emmeans(mod_success_false_MW_G_sub, ~Position+SMOOTH+Diversity+DC+Position:SMOOTH+Diversity:Position,
                     type = "response") %>%
               as_tibble(), METHOD = "FALSE")) %>%
  mutate(PD = paste(Position,Diversity, sep = " - ")) %>%
  mutate(PD = as.factor(PD)) %>%
  mutate(PD = fct_recode(PD,
                         "HR - R" = "high density level - high richness",
                         "LR - R" = "high density level - low richness",
                         "HR - L" = "low density level - high richness",
                         "LR - L" = "low density level - low richness",)) %>%
  ggplot(aes(y=prob,x=PD, color=METHOD,fill=METHOD, ymin=lower.CL, ymax=upper.CL))+
  facet_grid(DC~SMOOTH)+
  geom_hline(yintercept = seq(0,1,0.25), color="gray90",size=0.1)+
  geom_errorbar(width=0.2, color="gray60", position = position_dodge(width = 0.5),size=0.1)+
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
  #geom_point(shape=15,position = position_dodge(width = 0.5))+
  coord_cartesian(ylim = c(0,1))+
  scale_fill_manual("Position in sequence", labels=c("Focal area (correct detection)","Outside of focal area (false positive)"),
                    values = c("darkseagreen","coral"))+
  theme(axis.text.x = element_text(angle = -45, hjust = -0.05, vjust = 1))+
  labs(y="",
       x= ""),
nrow = 1, common.legend = T, legend = "top", labels = c("A","B"), widths = c(0.8,1)
) %>%
  annotate_figure(left = "Proportion of peak detection",
                  bottom = "Dataset type")

ggsave("METHOD_RESULTS//FIG2_v01.pdf",
       width = 20, height = 12, units = "cm")


emmeans(mod_success_focus, ~ PEAK,
        type = "response")


emmeans(mod_success_false, ~ PEAK,
        type = "response")


emmeans(mod_success_focus, ~ WU,
        type = "response")


emmeans(mod_success_false, ~ WU,
        type = "response")

emmeans(mod_success_focus, ~ WU*PEAK,
        type = "response") %>%
  as_tibble() %>%
  filter(PEAK == "PEAK.G")

emmeans(mod_success_false, ~ WU*PEAK,
        type = "response") %>%
  as_tibble() %>%
  filter(PEAK == "PEAK.G")

emmeans(mod_success_focus_MW_G_sub, ~ Position,
        type = "response")

emmeans(mod_success_false_MW_G_sub, ~Position,
        type = "response")


emmeans(mod_success_focus_MW_G_sub, ~ Diversity,
        type = "response")

emmeans(mod_success_false_MW_G_sub, ~Diversity,
        type = "response")

emmeans(mod_success_focus_MW_G_sub, ~ DC,
        type = "response")
emmeans(mod_success_false_MW_G_sub, ~DC,
        type = "response")

emmeans(mod_success_focus_MW_G_sub, ~ SMOOTH,
        type = "response")

##############
#   FIG 3
#############

get_dominant_pollen_taxa <- function(data, N=5){
  x = data$filtered.counts[-1]
  y = x/rowSums(x)*100
  z = colSums(y) %>% sort(decreasing =T)
  return(names(z)[1:N])
}

#SITE A


which(Hope_smooth$dataset.id %in%  17334 )

data_site_A <- list(dataset.id = Hope_smooth$dataset.id[[413]],
                    filtered.counts = Hope_smooth$filtered.counts[[413]],
                    list_ages = Hope_smooth$list_ages[[413]])

data_site_A_dom <- get_dominant_pollen_taxa(data_site_A)

data_site_A$filtered.counts %>%
  as_tibble() %>%
  dim()

#Site B

which(Hope_smooth$dataset.id %in%  4012 )


data_site_B <- list(dataset.id = Hope_smooth$dataset.id[[84]],
                    filtered.counts = Hope_smooth$filtered.counts[[84]],
                    list_ages = Hope_smooth$list_ages[[84]])

data_site_B_dom <- get_dominant_pollen_taxa(data_site_B)

data_site_B$filtered.counts %>%
  as_tibble() %>%
  dim()


# SITE c

which(Hope_smooth$dataset.id %in%  40951 )


data_site_C <- list(dataset.id = Hope_smooth$dataset.id[[273]],
                    filtered.counts = Hope_smooth$filtered.counts[[273]],
                    list_ages = Hope_smooth$list_ages[[273]])

data_site_C_dom <- get_dominant_pollen_taxa(data_site_C)

data_site_C$filtered.counts %>%
  as_tibble() %>%
  dim()


# Site D


which(Hope_smooth$dataset.id %in%  45314 )

data_site_D <- list(dataset.id = Hope_smooth$dataset.id[[299]],
                    filtered.counts = Hope_smooth$filtered.counts[[299]],
                    list_ages = Hope_smooth$list_ages[[299]])

data_site_D_dom <- get_dominant_pollen_taxa(data_site_D)

data_site_D$filtered.counts %>%
  as_tibble() %>%
  dim()


#Pollen taxa table

common_taxa<- c(data_site_A_dom,
                data_site_B_dom,
                data_site_C_dom,
                data_site_D_dom) %>%
  unique() %>%
  sub("/",".",.) %>%
  sub("-",".",.) %>%
  sub(")",".",.) %>%
  sub(".\\(","..",.) 

library (RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
Palette.1<- getPalette(length(common_taxa))
names(Palette.1)<- sort(common_taxa)  


fc_get_pollen_data <- function (data, sm.type, Common.list)
{
  # remove the sample ID
  if (is.numeric(unlist(data$filtered.counts[,1]))==F){
    data$filtered.counts <- data$filtered.counts[,-1]
  }
  
  data.ext <-  fc_extract_data(data$filtered.counts,
                               data$list_ages) %>%
    fc_smooth_pollen_data(.,sm.type = sm.type,
                          N.points = 5,
                          grim.N.max = 9,
                          range.age.max = 500) %>%
    fc_check_data(.,proportion = T)
  
  plot.data <- data.ext$Pollen %>%
    select(any_of(Common.list)) %>%
    rownames_to_column() %>%
    pivot_longer(cols = -c(rowname)) %>%
    rename(sample.id = rowname) %>%
    inner_join(.,data.ext$Age, by="sample.id")
  
  return (plot.data)
}


data_site_A_pollen <-fc_get_pollen_data(data_site_A, sm.type = "none",common_taxa)

data_site_B_pollen <-fc_get_pollen_data(data_site_B, sm.type = "none", common_taxa)

data_site_C_pollen <-fc_get_pollen_data(data_site_C, sm.type = "none",common_taxa)

data_site_D_pollen <-fc_get_pollen_data(data_site_D, sm.type = "none",common_taxa)


# Rate of Change

data_site_A_RoC_levels <- fc_R_ratepol(data.source.pollen = data_site_A$filtered.counts,
                                       data.source.age = data_site_A$list_ages,
                                       sm.type = "age.w", 
                                       N.points = 5,
                                       range.age.max = 500, 
                                       grim.N.max = 9,
                                       Working.Unit = "levels",
                                       BIN.size = 500,
                                       N.shifts = 5,
                                       rand = 10000,
                                       standardise = T, 
                                       S.value = 150, 
                                       DC = "chisq",
                                       interest.treshold = age_lim,
                                       Debug = F)

data_site_A_RoC_BINs <- fc_R_ratepol(data.source.pollen = data_site_A$filtered.counts,
                                     data.source.age = data_site_A$list_ages,
                                     sm.type = "age.w", 
                                     N.points = 5,
                                     range.age.max = 500, 
                                     grim.N.max = 9,
                                     Working.Unit = "BINs",
                                     BIN.size = 500,
                                     N.shifts = 5,
                                     rand = 10000,
                                     standardise = T, 
                                     S.value = 150, 
                                     DC = "chisq",
                                     interest.treshold = age_lim,
                                     Debug = F)

data_site_A_RoC_MW <- fc_R_ratepol(data.source.pollen = data_site_A$filtered.counts,
                                   data.source.age = data_site_A$list_ages,
                                   sm.type = "age.w", 
                                   N.points = 5,
                                   range.age.max = 500, 
                                   grim.N.max = 9,
                                   Working.Unit = "MW",
                                   BIN.size = 500,
                                   N.shifts = 5,
                                   rand = 10000,
                                   standardise = T, 
                                   S.value = 150, 
                                   DC = "chisq",
                                   interest.treshold = age_lim,
                                   Debug = F)

data_Site_B_RoC_levels <- fc_R_ratepol(data.source.pollen = data_site_B$filtered.counts,
                                       data.source.age = data_site_B$list_ages,
                                       sm.type = "shep", 
                                       N.points = 5,
                                       range.age.max = 500, 
                                       grim.N.max = 9,
                                       Working.Unit = "levels",
                                       BIN.size = 500,
                                       N.shifts = 5,
                                       rand = 10000,
                                       standardise = T, 
                                       S.value = 150, 
                                       DC = "chord",
                                       interest.treshold = age_lim,
                                       Debug = F)

data_Site_B_RoC_BINs <- fc_R_ratepol(data.source.pollen = data_site_B$filtered.counts,
                                     data.source.age = data_site_B$list_ages,
                                     sm.type = "shep", 
                                     N.points = 5,
                                     range.age.max = 500, 
                                     grim.N.max = 9,
                                     Working.Unit = "BINs",
                                     BIN.size = 500,
                                     N.shifts = 5,
                                     rand = 10000,
                                     standardise = T, 
                                     S.value = 150, 
                                     DC = "chord",
                                     interest.treshold = age_lim,
                                     Debug = F)

data_Site_B_RoC_MW <- fc_R_ratepol(data.source.pollen = data_site_B$filtered.counts,
                                   data.source.age = data_site_B$list_ages,
                                   sm.type = "shep", 
                                   N.points = 5,
                                   range.age.max = 500, 
                                   grim.N.max = 9,
                                   Working.Unit = "MW",
                                   BIN.size = 500,
                                   N.shifts = 5,
                                   rand = 10000,
                                   standardise = T, 
                                   S.value = 150, 
                                   DC = "chord",
                                   interest.treshold = age_lim,
                                   Debug = F)

data_Site_C_RoC_levels <- fc_R_ratepol(data.source.pollen = data_site_C$filtered.counts,
                                       data.source.age = data_site_C$list_ages,
                                       sm.type = "shep", 
                                       N.points = 5,
                                       range.age.max = 500, 
                                       grim.N.max = 9,
                                       Working.Unit = "levels",
                                       BIN.size = 500,
                                       N.shifts = 5,
                                       rand = 10000,
                                       standardise = T, 
                                       S.value = 150, 
                                       DC = "chord",
                                       interest.treshold = age_lim,
                                       Debug = F)

data_Site_C_RoC_BINs <- fc_R_ratepol(data.source.pollen = data_site_C$filtered.counts,
                                     data.source.age = data_site_C$list_ages,
                                     sm.type = "shep", 
                                     N.points = 5,
                                     range.age.max = 500, 
                                     grim.N.max = 9,
                                     Working.Unit = "BINs",
                                     BIN.size = 500,
                                     N.shifts = 5,
                                     rand = 10000,
                                     standardise = T, 
                                     S.value = 150, 
                                     DC = "chord",
                                     interest.treshold = age_lim,
                                     Debug = F)

data_Site_C_RoC_MW <- fc_R_ratepol(data.source.pollen = data_site_C$filtered.counts,
                                   data.source.age = data_site_C$list_ages,
                                   sm.type = "shep", 
                                   N.points = 5,
                                   range.age.max = 500, 
                                   grim.N.max = 9,
                                   Working.Unit = "MW",
                                   BIN.size = 500,
                                   N.shifts = 5,
                                   rand = 10000,
                                   standardise = T, 
                                   S.value = 150, 
                                   DC = "chord",
                                   interest.treshold = age_lim,
                                   Debug = F)

data_Site_D_RoC_levels <- fc_R_ratepol(data.source.pollen = data_site_D$filtered.counts,
                                       data.source.age = data_site_D$list_ages,
                                       sm.type = "shep", 
                                       N.points = 5,
                                       range.age.max = 500, 
                                       grim.N.max = 9,
                                       Working.Unit = "levels",
                                       BIN.size = 500,
                                       N.shifts = 5,
                                       rand = 10000,
                                       standardise = T, 
                                       S.value = 150, 
                                       DC = "chord",
                                       interest.treshold = age_lim,
                                       Debug = F)

data_Site_D_RoC_BINs <- fc_R_ratepol(data.source.pollen = data_site_D$filtered.counts,
                                     data.source.age = data_site_D$list_ages,
                                     sm.type = "shep", 
                                     N.points = 5,
                                     range.age.max = 500, 
                                     grim.N.max = 9,
                                     Working.Unit = "BINs",
                                     BIN.size = 500,
                                     N.shifts = 5,
                                     rand = 10000,
                                     standardise = T, 
                                     S.value = 150, 
                                     DC = "chord",
                                     interest.treshold = age_lim,
                                     Debug = F)


data_Site_D_RoC_MW <- fc_R_ratepol(data.source.pollen = data_site_D$filtered.counts,
                                   data.source.age = data_site_D$list_ages,
                                   sm.type = "shep", 
                                   N.points = 5,
                                   range.age.max = 500, 
                                   grim.N.max = 9,
                                   Working.Unit = "MW",
                                   BIN.size = 500,
                                   N.shifts = 5,
                                   rand = 10000,
                                   standardise = T, 
                                   S.value = 150, 
                                   DC = "chord",
                                   interest.treshold = age_lim,
                                   Debug = F)

# FIGURES 

FIG3_Site_A_1 <- data_site_A$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = c(0,3e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  geom_rug(sides = "b")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,3e-4))+
  labs(x="Age (cal yr BP)",
       y="Density of samples"
  )+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))


FIG3_Site_B_1 <- data_site_B$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = c(0,3e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  geom_rug(sides = "b")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,3e-4))+
  labs(x="Age (cal yr BP)",
       y="Density of samples"
  )+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))

FIG3_Site_C_1 <- data_site_C$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = c(0,3e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  geom_rug(sides = "b")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,3e-4))+
  labs(x="Age (cal yr BP)",
       y="Density of samples"
  )+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))


FIG3_Site_D_1 <- data_site_D$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = c(0,3e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  geom_rug(sides = "b")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  labs(x="Age (cal yr BP)",
       y="Density of samples"
  )+
  theme(
    #    axis.ticks.x = element_blank(),
    #    axis.text.x = element_blank(),
    #    axis.title.x = element_blank()
  )+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))


library(cowplot)
my_legend <- cowplot::get_legend(ggplot(data = data.frame(NAME=common_taxa,X=1), aes(x=X, fill=NAME))+
                                   geom_density(alpha=1/3)+
                                   theme_classic()+
                                   theme(legend.position = "bottom")+
                                   scale_fill_manual("pollen taxa",values = Palette.1))

plot(my_legend)


FIG3_Site_A_2 <-  data_site_A_pollen %>%
  bind_rows(.,data.frame(sample.id=NA,name=common_taxa,value=0,depth=NA, age = 10e3, newage=NA)) %>%
  ggplot(aes( y=value, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,1))+
  geom_hline(yintercept = c(0,1), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
              color="gray20", alpha=1/5, size=0.1)+
  scale_fill_manual("pollen taxa",values = Palette.1)+
  labs(
    x="Age (cal yr BP)",
    y="% of pollen grains"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(angle = 45)
        #        strip.text = element_blank()
  )+
  coord_flip(xlim=c(age_lim,0), ylim = c(0,1))+
  facet_wrap(~name, ncol=length(common_taxa))

FIG3_Site_B_2 <-  data_site_B_pollen %>%
  bind_rows(.,data.frame(sample.id=NA,name=common_taxa,value=0,depth=NA, age = 10e3, newage=NA)) %>%
  ggplot(aes( y=value, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,1))+
  geom_hline(yintercept = c(0,1), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
              color="gray20", alpha=1/5, size=0.1)+
  scale_fill_manual("pollen taxa",values = Palette.1, drop=FALSE)+
  labs(
    x="Age (cal yr BP)",
    y="% of pollen grains"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        strip.background = element_blank(),
        #        strip.text = element_text(angle = 45)
        strip.text = element_blank()
  )+
  coord_flip(xlim=c(age_lim,0), ylim = c(0,1))+
  facet_wrap(~name, ncol=length(common_taxa))

FIG3_Site_C_2 <-  data_site_C_pollen %>%
  bind_rows(.,data.frame(sample.id=NA,name=common_taxa,value=0,depth=NA, age = 10e3, newage=NA)) %>%
  ggplot(aes( y=value, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,1))+
  geom_hline(yintercept = c(0,1), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
              color="gray20", alpha=1/5, size=0.1)+
  scale_fill_manual("pollen taxa",values = Palette.1, drop=FALSE)+
  labs(
    x="Age (cal yr BP)",
    y="% of pollen grains"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        strip.background = element_blank(),
        #        strip.text = element_text(angle = 45)
        strip.text = element_blank()
  )+
  coord_flip(xlim=c(age_lim,0), ylim = c(0,1))+
  facet_wrap(~name, ncol=length(common_taxa))

FIG3_Site_D_2 <-  data_site_D_pollen %>%
  bind_rows(.,data.frame(sample.id=NA,name=common_taxa,value=0,depth=NA, age = 10e3, newage=NA)) %>%
  ggplot(aes( y=value, 
              x= age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,1))+
  geom_hline(yintercept = c(0,1), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=rep(0,length(value)),ymax=value, fill=name), 
              color="gray20", alpha=1/5, size=0.1)+
  scale_fill_manual("pollen taxa",values = Palette.1)+
  labs(
    x="Age (cal yr BP)",
    y="% of pollen grains"
  )+
  theme(legend.position = "none",
        #        axis.ticks.x = element_blank(),
        #        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        #        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        strip.background = element_blank(),
        #       strip.text = element_text(angle = 45)
        strip.text = element_blank()
  )+
  coord_flip(xlim=c(age_lim,0), ylim = c(0,1))+
  facet_wrap(~name, ncol=length(common_taxa))



FIG3_Site_A_3_levels <- data_site_A_RoC_levels %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  #geom_hline(yintercept = seq(from=0,to=5, by=2), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )



FIG3_Site_A_3_BINs <-  data_site_A_RoC_BINs %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  #geom_hline(yintercept = seq(from=0,to=2, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )
FIG3_Site_A_3_MW <-  data_site_A_RoC_MW %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  #geom_hline(yintercept = seq(from=0,to=2, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )

FIG3_Site_B_3_levels <-  data_Site_B_RoC_levels %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  # geom_hline(yintercept = seq(from=0,to=5, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )
FIG3_Site_B_3_BINs <-  data_Site_B_RoC_BINs %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  #geom_hline(yintercept = seq(from=0,to=2, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )

FIG3_Site_B_3_MW <-  data_Site_B_RoC_MW %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  # geom_hline(yintercept = seq(from=0,to=2, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )

FIG3_Site_C_3_levels <-data_Site_C_RoC_levels %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  # geom_hline(yintercept = seq(from=0,to=5, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )

FIG3_Site_C_3_BINs <-  data_Site_C_RoC_BINs %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  #geom_hline(yintercept = seq(from=0,to=2, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )

FIG3_Site_C_3_MW <-  data_Site_C_RoC_MW %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  #geom_hline(yintercept = seq(from=0,to=2, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )


FIG3_Site_D_3_levels <- data_Site_D_RoC_levels %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  # geom_hline(yintercept = seq(from=0,to=5, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        #        axis.ticks.x = element_blank(),
        #        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        #        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )

FIG3_Site_D_3_BINs <-  data_Site_D_RoC_BINs %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  # geom_hline(yintercept = seq(from=0,to=2, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        #        axis.ticks.x = element_blank(),
        #        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        #        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )
FIG3_Site_D_3_MW <-  data_Site_D_RoC_MW %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
  # geom_hline(yintercept = seq(from=0,to=2, by=0.5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BP)",
    y="Rate-of-Change score"
  )+
  theme(legend.position = "none",
        #        axis.ticks.x = element_blank(),
        #        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        #        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
  )


data_site_A_RoC_MW %>%
  filter(PEAK == T)

data_Site_B_RoC_MW %>%
  filter(PEAK == T)

data_Site_C_RoC_MW %>%
  filter(PEAK == T)

data_Site_D_RoC_MW %>%
  filter(PEAK == T)



library(cowplot)
FIG3_Site_A <- plot_grid(FIG3_Site_A_1, FIG3_Site_A_2, FIG3_Site_A_3_levels,FIG3_Site_A_3_BINs,FIG3_Site_A_3_MW,
                         align = "h",axis = "bt", ncol = 5, rel_widths = c(0.8,3,0.8,0.8,0.8))
FIG3_Site_B <- plot_grid(FIG3_Site_B_1, FIG3_Site_B_2, FIG3_Site_B_3_levels,FIG3_Site_B_3_BINs,FIG3_Site_B_3_MW,
                         align = "h",axis = "bt", ncol = 5, rel_widths = c(0.8,3,0.8,0.8,0.8))
FIG3_Site_C <- plot_grid(FIG3_Site_C_1, FIG3_Site_C_2, FIG3_Site_C_3_levels,FIG3_Site_C_3_BINs,FIG3_Site_C_3_MW,
                         align = "h",axis = "bt", ncol = 5, rel_widths = c(0.8,3,0.8,0.8,0.8))
FIG3_Site_D <- plot_grid(FIG3_Site_D_1, FIG3_Site_D_2, FIG3_Site_D_3_levels,FIG3_Site_D_3_BINs,FIG3_Site_D_3_MW,
                         align = "h",axis = "bt", ncol = 5, rel_widths = c(0.8,3,0.8,0.8,0.8))


FIG3_Site_comparison <- ggarrange(
  FIG3_Site_A,
  FIG3_Site_B,
  FIG3_Site_C,
  FIG3_Site_D,
  #labels = c(17334,"","",40951,"","",4012,"",""),
  labels = c("A","B","C","D"),
  heights = c(1.8,1,1,1.2),
  ncol = 1, nrow = 4, legend = "none")
FIG3_Site_comparison



ggsave("METHOD_RESULTS//FIG3_v02.pdf",
       width = 20, height = 22, units = "cm")





############
#   FIG 4
############

fc_calculate_RoC_comparison <- function(data, Working.Unit, BIN.size, N.shifts, rand ,peek, interest.treshold){
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  
  for(i in 1:20)
  {
    data.temp<- fc_R_ratepol( data.source.pollen =  data$filtered.counts,
                              data.source.age = data$list_ages,
                              sm.type = performance.smooth[i],
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              Working.Unit = Working.Unit,
                              BIN.size = BIN.size,
                              N.shifts = N.shifts,
                              rand = rand,
                              standardise = F, 
                              S.value = 150 ,
                              DC = performance.DC[i],
                              interest.treshold = interest.treshold,
                              Peak = peek,
                              Debug = F) %>%
      as_tibble()
    
    if(peek == "all"){
      # PEAK detection 
      # Median peak treshold
      # treshold for RoC peaks is set as median of all RoC in dataset
      r.treshold <- median(data.temp$ROC)
      # mark peaks which have 95% quantile above the treshold asPeak.treshold
      data.temp$PEAK.T <- data.temp$ROC.dw > r.treshold
      
      # GAM  
      # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
      pred.gam <-  predict.gam(gam(ROC~s(AGE,k=3), data = data.temp, family = "Gamma",
                                   correlation = corCAR1(form = ~ AGE), method = "REML"), type="response")
      pred.gam.diff <- data.temp$ROC - pred.gam
      data.temp$PEAK.G <- (pred.gam.diff) > 1.5*sd(pred.gam.diff)
      
      # SNI  
      # set moving window of 5 times higher than average distance between samples
      mean.age.window <- 5 * mean( diff(data.temp$AGE) )
      # calculate SNI (singal to noise ratio)
      SNI.calc <- CharSNI(data.frame(data.temp$AGE, data.temp$ROC, pred.gam),mean.age.window)
      # mark points with SNI higher than 3
      data.temp$PEAK.S <- SNI.calc$SNI > 3 & data.temp$ROC > pred.gam
    }
    
    data.temp.sum<-data.temp %>%
      mutate(SMOOTH = performance.smooth[i],
             DC = performance.DC[i]) %>%
      select(SMOOTH, DC, Working_Unit ,AGE, ROC,ROC.up, ROC.dw, PEAK)
    
    if(i == 1){
      res.tibble <- data.temp.sum} else {
        res.tibble <- rbind(res.tibble,data.temp.sum)
      }
  }
  return(res.tibble)
}


data_example <- list(dataset.id = tibble_Europe2$dataset.id[[2]],
                     filtered.counts = tibble_Europe2$filtered.counts[[2]],
                     list_ages = tibble_Europe2$list_ages[[2]])

data_example_MW <- fc_calculate_RoC_comparison(data_example,
                                               Working.Unit = "MW",
                                               BIN.size = 500,
                                               N.shifts = 5,
                                               rand = 10000,
                                               peek = "GAM",
                                               interest.treshold =  age_lim)
data_example_MW.fresh <- data_example_MW
data_example_MW <- data_example_MW.fresh


data_example_MW <- within(data_example_MW, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
levels(data_example_MW$DC) <- c("Euc","Euc.sd","Chord","Chisq")

data_example_MW <- within(data_example_MW, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))
levels(data_example_MW$SMOOTH)<-c("None","M.avg","Grimm","Age.w","Shep")

data_example_MW[data_example_MW$DC == "Chisq" & data_example_MW$SMOOTH == "Age.w",]$PEAK[1] <- F


FIG1_visual_example_MW <- data_example_MW %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0))+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  xlab("Age (cal yr BP)")+ylab("Rate of Change score")+
  facet_grid(SMOOTH~DC, scales = "free_x")

FIG1_visual_example_MW

ggsave("METHOD_RESULTS/FIG4_v01.pdf",
       plot = FIG1_visual_example_MW,
       width = 16, height = 12, units = "cm")


data_example_BIN <- fc_calculate_RoC_comparison(data_example,
                                                Working.Unit = "BINs",
                                                BIN.size = 500,
                                                rand = 10000,
                                                interest.treshold =  age_lim)
data_example_BIN <- within(data_example_BIN, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
levels(data_example_BIN$DC) <- c("Euc","Euc.sd","Chord","Chisq")
data_example_BIN <- within(data_example_BIN, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))
levels(data_example_BIN$SMOOTH)<-c("None","M.avg","Grimm","Age.w","Shep")

FIG1_visual_example_BIN <- data_example_BIN %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0))+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK.T==T),color="blue", size=2, shape=1, alpha=2/3)+
  geom_point(data = . %>% filter(PEAK.G==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_point(data = . %>% filter(PEAK.S==T ),color="red", size=2, shape=8, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  xlab("Age (cal yr BP)")+ylab("Rate of Change score")+
  facet_grid(SMOOTH~DC, scales = "free_x")

FIG1_visual_example_BIN

ggsave("~/RESULTS/Methods/FIN/FIG1_visual_example_BIN.pdf",
       plot = FIG1_visual_example_BIN,
       width = 20, height = 12, units = "cm")



###############
# SUPLEMENTARY 
###############

##########
# FIG S1 #
########## 

# Example of simulation of enviroemntal data

time= tibble_Europe2$list_ages[[2]]$ages$age
nforc=4;
mean=100; 
sdev=.15; 
nprox=50; 
var=20;
range=15;
manual.edit = T;
breaks=c(2000,3000);
breaks=c(5500,6500);
jitter = T;
rarity=T;

forcing <- array(0, dim=c(length(time), nforc))
forcing[1,] <- rnorm(nforc, mean, sdev)
for(i in 2:length(time))
  forcing[i,] <- rnorm(nforc, forcing[i-1,], sdev)
for(l in 1:( length(breaks)-1 ) )
{
  if(l%%2 == 1) # odd
  {
    forcing[time>breaks[l] & time<breaks[l+1],] <-  forcing[time>breaks[l] & time<breaks[l+1],] * (1+sdev)  
  }
  
  if(l%%2 == 0) # even
  {
    forcing[time>breaks[l] & time<breaks[l+1],] <-  forcing[time>breaks[l] & time<breaks[l+1],] * (1-sdev) 
  }
}

#forcing <- saved.f
#saved.f <-forcing 

forcing<- apply(forcing,2, FUN = function(x) {
  low <- lowess(x,f=.05,iter=100)
  return(low$y)
}) 

# choose random optima and ranges for the biota
ecology <- c()
for(i in 1:nprox)
  ecology[[i]] <- list(mean=rnorm(nforc, mean, var), sd=rgamma(nforc, range, 1))

# reactions of the biota to the environmental changes
proxies <- array(1, dim=c(length(time), nprox))
o <- c()

for(i in 1:nprox)
{
  for(j in 1:nforc)
    proxies[,i] <- proxies[,i] * dnorm(forcing[,j], ecology[[i]]$mean[j], ecology[[i]]$sd[j])
  #o[i] <- weighted.mean(time, proxies[,i])
}

# order taxa by abundance
o <- order(colSums(proxies), decreasing=TRUE)
proxies <- proxies[,o]


# decrease the abundances of rare taxa
if(rarity==T)
{
  for(i in 1:ncol(proxies)) {
    proxies[,i]<- (proxies[,i] / max(1,runif(1, min = i-1, max=i)) )
  }
}


# jitter the resul the pollen data
if(jitter==T)
{
  proxies<- apply(proxies,2,FUN= function(x) jitter(x,factor = 1.5, amount = 0))
  proxies[proxies < 0] <- 0
}

# return 
data.source.age<- list(ages=data.frame(sample.id =c(1:length(time)), age=time),
                       age_position= matrix(time, nrow = 1))
data.source.pollen <- as.data.frame(proxies)
row.names(data.source.pollen) <- c(1:length(time))



Supplementary_F1a <-ggarrange(as.data.frame(forcing) %>%
                               mutate(AGE = time) %>%
                               pivot_longer(., cols = -c(AGE)) %>%
                               arrange(AGE,value) %>%
                               ggplot(aes(x=AGE, y= value))+
                               geom_vline(xintercept = breaks, color="gray80", size=0.1)+
                               geom_line(aes(color=name))+
                               theme_classic()+
                               coord_flip(xlim=c(8000,0))+
                               scale_x_continuous(trans = "reverse")+
                               theme(
                                 axis.ticks.y = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.title.y = element_blank(),
                                     #axis.ticks.x = element_blank(),
                                     #axis.text.x = element_blank(),
                                     legend.position = "none"
                                     )+
                                #ylab("")+
                                ylab("Value of env. variable")+
                                xlab("Age (cal yr BP)"),
                             fc_extract_data(data.source.pollen, data.source.age) %>%
                               fc_smooth_pollen_data("none") %>%
                               fc_check_data(., proportion = T) %>%
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
                               xlab("")+
                               #ylab("")+
                               ylab("Proportion of pollen grains")+
                               theme(
                                     axis.ticks.y = element_blank(),
                                     axis.text.y = element_blank(),
                                     #axis.ticks.x = element_blank(),
                                     #axis.text.x = element_blank(),
                                     legend.position = "none"),
                             ncol=2, align = "h")


Supplementary_F1a


Supplementary_F1 <- ggarrange(
  Supplementary_F1a,Supplementary_F1b,
  nrow = 2, labels=c("A","B")
)


Supplementary_F1_fin <- annotate_figure(Supplementary_F1, left= "Age (cal yr BP)")

Supplementary_F1_fin

ggsave("~/RESULTS/Methods/FIN/Supplementary_F1.pdf",
       plot = Supplementary_F1_fin,
       height = 10, width = 15, units="cm")

##########
# FIG S2 #
########## 
ggarrange(ggarrange(emmeans(mod_success_focus_MW_G_sub, ~ DC,
                            type = "response") %>%
                      as_tibble() %>%
                      ggplot(aes(y=prob, x=DC,color=DC,fill=DC, ymin=lower.CL, ymax= upper.CL))+
                      geom_hline(yintercept = seq(0,1,0.1), color="gray90",size=0.1)+
                      geom_errorbar(width=0.2, color="gray60",position = position_dodge(width = 0.5),size=0.1)+
                      geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
                      #geom_point(shape=15,position = position_dodge(width = 0.5))+
                      coord_cartesian(ylim = c(0.7,1))+
                      scale_fill_manual(values = Color.legen_DC)+
                      scale_color_manual(values = Color.legen_DC)+
                      labs(y="",x=""),
                    emmeans(mod_success_false_MW_G_sub, ~ DC,
                            type = "response") %>%
                      as_tibble() %>%
                      ggplot(aes(y=prob, x=DC,color=DC,fill=DC, ymin=lower.CL, ymax= upper.CL))+
                      geom_hline(yintercept = seq(0,1,0.025), color="gray90",size=0.1)+
                      geom_errorbar(width=0.2, color="gray60",position = position_dodge(width = 0.5),size=0.1)+
                      geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
                      #geom_point(shape=15,position = position_dodge(width = 0.5))+
                      coord_cartesian(ylim = c(0,0.1))+
                      scale_fill_manual(values = Color.legen_DC)+
                      scale_color_manual(values = Color.legen_DC)+
                      labs(y="",x=""),
                    nrow = 1, common.legend = T, legend = "none",
                    labels = c("Correct detections","False positives")
),
ggarrange(emmeans(mod_success_focus_MW_G_sub, ~ SMOOTH,
                  type = "response") %>%
            as_tibble() %>%
            ggplot(aes(y=prob, x=SMOOTH,color=SMOOTH,fill=SMOOTH, ymin=lower.CL, ymax= upper.CL))+
            geom_hline(yintercept = seq(0,1,0.1), color="gray90",size=0.1)+
            geom_errorbar(width=0.2, color="gray60", position = position_dodge(width = 0.5),size=0.1)+
            geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
            #geom_point(shape=15,position = position_dodge(width = 0.5))+
            coord_cartesian(ylim = c(0.7,1))+
            scale_color_manual(values = Color.legen_smooth)+
            scale_fill_manual(values = Color.legen_smooth)+
            labs(y="",x=""),
          emmeans(mod_success_false_MW_G_sub, ~ SMOOTH,
                  type = "response") %>%
            as_tibble() %>%
            ggplot(aes(y=prob, x=SMOOTH,color=SMOOTH,fill=SMOOTH, ymin=lower.CL, ymax= upper.CL))+
            geom_hline(yintercept = seq(0,1,0.05), color="gray90",size=0.1)+
            geom_errorbar(width=0.2, color="gray60", position = position_dodge(width = 0.5),size=0.1)+
            geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
            #geom_point(shape=15,position = position_dodge(width = 0.5))+
            coord_cartesian(ylim = c(0,0.2))+
            scale_color_manual(values = Color.legen_smooth)+
            scale_fill_manual(values = Color.legen_smooth)+
            labs(y="",x=""),
          nrow = 1, common.legend = T, legend = "none"),
ggarrange(emmeans(mod_success_focus_MW_G_sub, ~ Position ,
                  type = "response") %>%
            as_tibble() %>%
            ggplot(aes(y=prob, x=Position ,color=Position ,fill=Position , ymin=lower.CL, ymax= upper.CL))+
            geom_hline(yintercept = seq(0,1,0.1), color="gray90",size=0.1)+
            geom_errorbar(width=0.2, color="gray60", position = position_dodge(width = 0.5),size=0.1)+
            geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
            #geom_point(shape=15,position = position_dodge(width = 0.5))+
            coord_cartesian(ylim = c(0.5,1))+
            scale_color_manual(values = Color.legen_Position )+
            scale_fill_manual(values = Color.legen_Position )+
            labs(y="",x=""),
          emmeans(mod_success_false_MW_G_sub, ~ Position ,
                  type = "response") %>%
            as_tibble() %>%
            ggplot(aes(y=prob, x=Position ,color=Position ,fill=Position , ymin=lower.CL, ymax= upper.CL))+
            geom_hline(yintercept = seq(0,1,0.05), color="gray90",size=0.1)+
            geom_errorbar(width=0.2, color="gray60", position = position_dodge(width = 0.5),size=0.1)+
            geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
            #geom_point(shape=15,position = position_dodge(width = 0.5))+
            coord_cartesian(ylim = c(0,0.2))+
            scale_color_manual(values = Color.legen_Position )+
            scale_fill_manual(values = Color.legen_Position )+
            labs(y="",x=""),
          nrow = 1, common.legend = T, legend = "none"),
ggarrange(emmeans(mod_success_focus_MW_G_sub, ~ Diversity ,
                  type = "response") %>%
            as_tibble() %>%
            ggplot(aes(y=prob, x=Diversity ,color=Diversity ,fill=Diversity , ymin=lower.CL, ymax= upper.CL))+
            geom_hline(yintercept = seq(0,1,0.1), color="gray90",size=0.1)+
            geom_errorbar(width=0.2, color="gray60", position = position_dodge(width = 0.5),size=0.1)+
            geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
            #geom_point(shape=15,Diversity = position_dodge(width = 0.5))+
            coord_cartesian(ylim = c(0.7,1))+
            scale_color_manual(values = Color.legen_Diversity )+
            scale_fill_manual(values = Color.legen_Diversity )+
            labs(y="",x=""),
          emmeans(mod_success_false_MW_G_sub, ~ Diversity ,
                  type = "response") %>%
            as_tibble() %>%
            ggplot(aes(y=prob, x=Diversity ,color=Diversity ,fill=Diversity , ymin=lower.CL, ymax= upper.CL))+
            geom_hline(yintercept = seq(0,1,0.1), color="gray90",size=0.1)+
            geom_errorbar(width=0.2, color="gray60", position = position_dodge(width = 0.5),size=0.1)+
            geom_bar(stat = "identity", position = position_dodge(width = 0.5), color="gray60", orientation="x", width = 0.5,size=0.1)+
            #geom_point(shape=15,Diversity = position_dodge(width = 0.5))+
            coord_cartesian(ylim = c(0,0.4))+
            scale_color_manual(values = Color.legen_Diversity )+
            scale_fill_manual(values = Color.legen_Diversity )+
            labs(y="",x=""),
          nrow = 1, common.legend = T, legend = "none") ,
nrow = 4) %>% annotate_figure(left = "Proportion of peak detection")




ggsave("METHOD_RESULTS//FIG_S2_v01.pdf",
       width = 20, height = 20, units = "cm")



citation(package = "emmeans")




##########
# FIG S3 #
########## 


formula(glm.fin)


Success_supp_A <- ggarrange(
  data_success_sum %>%
    ungroup() %>%
    filter(PEAK == "PEAK.G") %>%
    filter(SEGMENT == "focus") %>%
    group_by(Position) %>%
    summarise(
      VALUE.M = mean(VALUE),
      SD = sd(VALUE),
      SE = SD/sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(x=Position, y=VALUE.M))+
    geom_bar(aes(fill=Position), stat = "identity", color="gray50")+
    geom_errorbar(aes(ymin=VALUE.M-SE, ymax=VALUE.M+SE), width=0.2, size=0.5, color="gray50")+
    theme_classic()+
    scale_fill_manual(values = Color.legen_Position)+
    coord_cartesian(ylim=c(0.4,1))+
    labs(x="Density of levels",
         y="Percentage of peak detection"
    )+
    theme(legend.position = "none", 
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  ,
  data_success_sum %>%
    ungroup() %>%
    filter(PEAK == "PEAK.G") %>%
    filter(SEGMENT == "focus") %>%
    group_by(SMOOTH) %>%
    summarise(
      VALUE.M = mean(VALUE),
      SD = sd(VALUE),
      SE = SD/sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(x=SMOOTH, y=VALUE.M))+
    geom_bar(aes(fill = SMOOTH),stat = "identity", color="gray50")+
    geom_errorbar(aes(ymin=VALUE.M-SE, ymax=VALUE.M+SE), width=0.2, size=0.5, color="gray50")+
    theme_classic()+
    scale_fill_manual(values = Color.legen_smooth)+
    coord_cartesian(ylim=c(0.4,1))+
    labs(x="Smoothing",
         y="Percentage of peak detection")+
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  ,
  data_success_sum %>%
    ungroup() %>%
    filter(PEAK == "PEAK.G") %>%
    filter(SEGMENT == "focus") %>%
    group_by(Diversity) %>%
    summarise(
      VALUE.M = mean(VALUE),
      SD = sd(VALUE),
      SE = SD/sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(x=Diversity, y=VALUE.M))+
    geom_bar(aes(fill=Diversity), stat = "identity", color="gray50")+
    geom_errorbar(aes(ymin=VALUE.M-SE, ymax=VALUE.M+SE), width=0.2, size=0.5, color="gray50")+
    theme_classic()+
    scale_fill_manual(values=Color.legen_Diversity)+
    coord_cartesian(ylim=c(0.4,1))+
    labs(x= "Richness of pollen taxa",
         y= "Percentage of peak detection")+
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  ,
  data_success_sum %>%
    ungroup() %>%
    filter(PEAK == "PEAK.G") %>%
    filter(SEGMENT == "focus") %>%
    group_by(DC) %>%
    summarise(
      VALUE.M = mean(VALUE),
      SD = sd(VALUE),
      SE = SD/sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(x=DC, y=VALUE.M))+
    geom_bar(aes(fill=DC),stat = "identity", color="gray50", )+
    geom_errorbar(aes(ymin=VALUE.M-SE, ymax=VALUE.M+SE), width=0.2, size=0.5, color="gray50")+
    theme_classic()+
    scale_fill_manual(values = Color.legen_DC)+
    coord_cartesian(ylim=c(0.4,1))+
    labs(x= "Dissimilarity coeficient",
         y= "Percentage of peak detection")+
    theme(legend.position = "none", 
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  , nrow = 1, ncol = 4
)

Success_supp_A

Success_supp_A_a <- annotate_figure(Success_supp_A, top = "Focal area (correct detection)")

Success_supp_A_a

formula(glm.fin.e)

Success_supp_B <- ggarrange(
  data_success_sum %>%
    ungroup() %>%
    filter(PEAK == "PEAK.G") %>%
    filter(SEGMENT == "empty") %>%
    group_by(Position) %>%
    summarise(
      VALUE.M = mean(VALUE),
      SD = sd(VALUE),
      SE = SD/sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(x=Position, y=VALUE.M))+
    geom_bar(aes(fill=Position), stat = "identity", color="gray50")+
    geom_errorbar(aes(ymin=VALUE.M-SE, ymax=VALUE.M+SE), width=0.2, size=0.5, color="gray50")+
    theme_classic()+
    scale_fill_manual(values = Color.legen_Position)+
    coord_cartesian(ylim=c(0,0.25))+
    labs(x="Density of levels",
         y="Percentage of peak detection"
    )+
    theme(legend.position = "none", 
          axis.title.y = element_blank())
  ,
  data_success_sum %>%
    ungroup() %>%
    filter(PEAK == "PEAK.G") %>%
    filter(SEGMENT == "empty") %>%
    group_by(SMOOTH) %>%
    summarise(
      VALUE.M = mean(VALUE),
      SD = sd(VALUE),
      SE = SD/sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(x=SMOOTH, y=VALUE.M))+
    geom_bar(aes(fill = SMOOTH),stat = "identity", color="gray50")+
    geom_errorbar(aes(ymin=VALUE.M-SE, ymax=VALUE.M+SE), width=0.2, size=0.5, color="gray50")+
    theme_classic()+
    scale_fill_manual(values = Color.legen_smooth)+
    coord_cartesian(ylim=c(0,0.25))+
    labs(x="Smoothing",
         y="Percentage of peak detection")+
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  ,
  data_success_sum %>%
    ungroup() %>%
    filter(PEAK == "PEAK.G") %>%
    filter(SEGMENT == "empty") %>%
    group_by(Diversity) %>%
    summarise(
      VALUE.M = mean(VALUE),
      SD = sd(VALUE),
      SE = SD/sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(x=Diversity, y=VALUE.M))+
    geom_bar(aes(fill=Diversity), stat = "identity", color="gray50")+
    geom_errorbar(aes(ymin=VALUE.M-SE, ymax=VALUE.M+SE), width=0.2, size=0.5, color="gray50")+
    theme_classic()+
    scale_fill_manual(values=Color.legen_Diversity)+
    coord_cartesian(ylim=c(0,0.25))+
    labs(x= "Richness of pollen taxa",
         y= "Percentage of peak detection")+
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  ,
  data_success_sum %>%
    ungroup() %>%
    filter(PEAK == "PEAK.G") %>%
    filter(SEGMENT == "empty") %>%
    group_by(DC) %>%
    summarise(
      VALUE.M = mean(VALUE),
      SD = sd(VALUE),
      SE = SD/sqrt(n())
    ) %>%
    ungroup() %>%
    ggplot(aes(x=DC, y=VALUE.M))+
    geom_bar(aes(fill=DC),stat = "identity", color="gray50", )+
    geom_errorbar(aes(ymin=VALUE.M-SE, ymax=VALUE.M+SE), width=0.2, size=0.5, color="gray50")+
    theme_classic()+
    scale_fill_manual(values = Color.legen_DC)+
    coord_cartesian(ylim=c(0,0.25))+
    labs(x= "Dissimilarity coeficient",
         y= "Percentage of peak detection")+
    theme(legend.position = "none", 
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  , nrow = 1, ncol = 4
)


Success_supp_B

Success_supp_B_a <- annotate_figure(Success_supp_B, top="Outside of focal area (false positive)")

Success_supp_B_a

Success_supp <- ggarrange(
  Success_supp_A_a,
  Success_supp_B_a,
  nrow = 2
)

Success_supp

Supplementary_F3 <- annotate_figure(Success_supp, left = "Proportion of peak detection")

Supplementary_F3

ggsave("~/RESULTS/Methods/FIN/Supplementary_F3.pdf",
       plot = Supplementary_F3,
       height = 12, width = 25, units="cm")


##########
# FIG S4 #
########## 



Supplementary_F4 <- data_success_sum %>%
  ungroup() %>%
  filter(PEAK == "PEAK.G") %>%
  filter(SEGMENT == "empty") %>%
  group_by(DC, SMOOTH) %>%
  summarise(
    VALUE.M = mean(VALUE),
    SD = sd(VALUE),
    SE = SD/sqrt(n())
  ) %>%
  ungroup() %>% 
  ggplot(aes(x=SMOOTH, y= VALUE.M, group= DC, fill= DC))+
  geom_bar(stat = "identity", position = position_dodge(), color="gray30")+
  geom_errorbar(aes(ymin = VALUE.M-SE, ymax= VALUE.M+SE), position = position_dodge(width=0.9), color="gray50",width=0.2, size=0.5)+
  #coord_cartesian(ylim=c(0.62, 0.75))+
  theme_classic()+
  scale_fill_manual(values = Color.legen_DC)+
  labs(x= "Smoothing",
       y= "Proportion of peak detection",
       title = "Outside of focal area (false positive)",
       fill="Disimilarity coeficient")+
  theme()

Supplementary_F4

data_success_sum %>%
  ungroup() %>%
  filter(PEAK == "PEAK.G") %>%
  filter(SEGMENT == "empty") %>%
  group_by(DC, SMOOTH) %>%
  summarise(
    VALUE.M = mean(VALUE),
    SD = sd(VALUE),
    SE = SD/sqrt(n())
  ) %>%
  ungroup() %>%
  View()


ggsave("~/RESULTS/Methods/FIN/Supplementary_F4.pdf",
       plot = Supplementary_F4,
       height = 12, width = 25, units="cm")


##########
# FIG S5 #
########## 

Supplementary_F5<- tibble_Europe2[c(2,45,224,50),] %>%
  ggplot(aes(y=lat, x=long, label = c("A","B","C","D")))+
  borders(fill = "gray90", colour = "gray60") +
  geom_point(size=3)+
  geom_text(nudge_x = 1,nudge_y = -1)+
  coord_fixed(xlim = c(-25,25), ylim=c(35,70))+
  theme_classic()+
  labs(x="longitude",
       y="lattitude")

Supplementary_F5

ggsave("~/RESULTS/Methods/FIN/Supplementary_F5.pdf",
       plot = Supplementary_F5,
       height = 10, width = 15, units="cm")

##########
# FIG S6 #
########## 

data_scheme_00 <- data_site_A$list_ages$ages %>%
  filter(age < 8000)


data_scheme_1_a <-fc_get_pollen_data(data_site_A, sm.type = "none",N.taxa = 1)
data_scheme_1_a <-data_scheme_1_a %>%
  filter(age < 8000)
data_scheme_1_b <-fc_get_pollen_data(data_site_A, sm.type = "grim",N.taxa = 1)
data_scheme_1_b <-data_scheme_1_b %>%
  filter(age < 8000)

dia01 <- ggplot()+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_cartesian (ylim = c(0.1,1), xlim = c(8e3,0))+
  labs(x="Age (yr BP)",
       title = "Smoothing of pollen data")+
  theme(axis.title.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.y = element_blank())+
  geom_line(data = data_scheme_1_a,aes(x=age, y=value), color="gray50")+
  geom_line(data = data_scheme_1_b,aes(x=age, y=value+0.5), color="gray30")+
  geom_segment(aes(x=4e3, xend=4e3, y=0.3, yend=0.7), arrow = arrow(length = unit(0.3, "cm")), color="gray30")


dia01

BINs_a <-seq(from=0,
             to=ceiling(max(data_scheme_01$age)),
             by=500)

BINs_b <-seq(from=333,
             to=ceiling(max(data_scheme_01$age)),
             by=500)

BINs_c <-seq(from=666,
             to=ceiling(max(data_scheme_01$age)),
             by=500)

dia02 <- ggplot()+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(x= "Age (yr BP",
       title = "Creation of time bins")+
  geom_segment(aes(x = BINs_a, xend = BINs_a, y = 1, yend = 2), color = "blue")+
  geom_segment(aes(x = BINs_a[-length(BINs_a)]+100, xend = BINs_a[-1]-100, y = 1.5, yend = 1.5),
               color = "blue")+
  geom_segment(aes(x = BINs_b, xend = BINs_b, y = 2.5, yend = 3.5), color = "green")+
  geom_segment(aes(x = BINs_b[-length(BINs_b)]+100, xend = BINs_b[-1]-100, y = 3, yend = 3),
               color = "green")+
  geom_segment(aes(x = BINs_c, xend = BINs_c, y = 4, yend = 5), color = "red")+
  geom_segment(aes(x = BINs_c[-length(BINs_c)]+100, xend = BINs_c[-1]-100, y = 4.5, yend = 4.5),
               color = "red")

dia02

select_bins_help_fc <- function(y){
  
  x <- data.frame(matrix(ncol = ncol(data_scheme_01), nrow = length(y)))
  names(x) <- names(data_scheme_01)
  
  for( i in 1:length(y)){
    
    selected.BIN <- y[i] #select teh bin
    
    subset.w <- data_scheme_01 %>%
      filter(age < selected.BIN+500 & age > selected.BIN)
    
    if (nrow(subset.w)>0) # If selected subset has at least one sample
    {
      
      subset.w$diff <- abs(subset.w$age-selected.BIN)
      suppressWarnings(x[i,] <- subset.w[subset.w$diff==min(subset.w$diff),c(1:4)] )
      
    }
  }
  
  return(x)
} 



data_scheme_02_a <- na.omit(select_bins_help_fc(y=BINs_a))
data_scheme_02_b <- na.omit(select_bins_help_fc(y=BINs_b))
data_scheme_02_c <- na.omit(select_bins_help_fc(y=BINs_c))

data_scheme_03 <- tibble(age=c((data_scheme_02_a$age[-length(data_scheme_02_a$age)]+data_scheme_02_a$age[-1])/2,
                               (data_scheme_02_b$age[-length(data_scheme_02_b$age)]+data_scheme_02_b$age[-1])/2,
                               (data_scheme_02_c$age[-length(data_scheme_02_c$age)]+data_scheme_02_c$age[-1])/2),
                         color = c(rep("blue",nrow(data_scheme_02_a)-1),
                                   rep("green",nrow(data_scheme_02_b)-1),
                                   rep("red",nrow(data_scheme_02_c)-1)),
                         VALUE = NA) %>%
  arrange(age)



for(i in 1:nrow(data_scheme_03)){
  
  if(i ==1){
    data_scheme_03$VALUE[i] <- 1
  } else {
    data_scheme_03$VALUE[i] <- data_scheme_03$VALUE[i-1] + runif(1,min = -0.1, max = 0.1)
    
  }
}


data_scheme_03$VALUE <- lowess(data_scheme_03$age,data_scheme_03$VALUE,f=.2,iter=100)$y

data_scheme_03_a<- data_scheme_03 %>%
  filter(color == "blue")

data_scheme_03_b<- data_scheme_03 %>%
  filter(color == "green") 

data_scheme_03_c<- data_scheme_03 %>%
  filter(color == "red") 


ggplot()+
  theme_classic()+
  theme (axis.text.y = element_blank(),
         axis.ticks.y = element_blank())+
  labs(x="Age (cal yr BP)",
       y= "Rate-of-Change score",
       title = "Summarisation results from all Mowing Windows")+
  coord_cartesian(xlim = c(0,8000))+
  geom_point(aes(x=data_scheme_03_a$age, y=data_scheme_03_a$VALUE), shape=18, color="blue", size=3)+
  geom_point(aes(x=data_scheme_03_b$age, y=data_scheme_03_b$VALUE), shape=18, color="green", size=3)+
  geom_point(aes(x=data_scheme_03_c$age, y=data_scheme_03_c$VALUE), shape=18, color="red", size=3)+
  geom_line(data=data_scheme_03, aes(x=age, y=VALUE), color="gray30", lty=2)


dia03_A1<- ggplot()+
  theme_classic()+
  theme (axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.title.x= element_blank())+
  labs(x="Age (cal yr BP)")+
  coord_cartesian(xlim = c(0,8000))+
  geom_point(data = data_scheme_01, aes(x=age, y=1.3), shape = 0 , size=2, color="gray20")+
  geom_segment(aes(x = 4e3, y = 1.4, xend = 4e3, yend = 1.7), arrow = arrow(length = unit(0.3, "cm")), color = "gray80")+
  geom_segment(aes(x = BINs_a, y = 1.75, xend = BINs_a, yend = 2.25), color = "blue")+
  geom_point(data = data_scheme_01, aes(x=age, y=2), shape = 0 , size=2, color="gray80")+
  geom_point(data= data_scheme_02_a, aes(x=age, y= 2), shape = 15, size= 2, color="gray20")+
  geom_segment(data= data_scheme_02_a, aes(x = age, y = 2.1, xend = age, yend = 2.9), color = "gray80", lty = 3)+
  geom_point(data= data_scheme_02_a, aes(x=age, y= 3), shape = 15, size= 2, color="blue")

dia03_A1

dia03_A2 <- ggplot()+
  theme_classic()+
  theme (
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  )+
  labs(x="Age (cal yr BP)")+
  coord_cartesian(xlim = c(0,8000))+
  geom_point(data= data_scheme_02_a, aes(x=age, y= 1), shape = 15, size= 2, color="blue")+
  geom_segment(aes(x = data_scheme_02_a$age[-length(data_scheme_02_a$age)]+100, 
                   xend = data_scheme_02_a$age[-1]-100,
                   y = 1.1, 
                   yend = 1.1), color = "gray50", size= 3)+
  geom_segment(aes(x = data_scheme_03_a$age, y = 1.4,
                   xend = data_scheme_03_a$age, yend = data_scheme_03_a$VALUE+0.9), color = "gray50", lty = 3)+
  #geom_line(aes(x=data_scheme_03_a$age, y=data_scheme_03_a$VALUE+1), lty=2)+
  geom_point(aes(x=data_scheme_03_a$age, y=data_scheme_03_a$VALUE+1), shape=18, color="blue", size=3)

dia03_A2

dia03_B1<-ggplot()+
  theme_classic()+
  theme (axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.title.x= element_blank())+
  labs(x="Age (cal yr BP)")+
  coord_cartesian(xlim = c(0,8000))+
  geom_point(data = data_scheme_01, aes(x=age, y=1.3), shape = 0 , size=2, color="gray20")+
  geom_segment(aes(x = 4e3, y = 1.4, xend = 4e3, yend = 1.7), arrow = arrow(length = unit(0.3, "cm")), color = "gray80")+
  geom_segment(aes(x = BINs_b, y = 1.75, xend = BINs_b, yend = 2.25), color = "green")+
  geom_point(data = data_scheme_01, aes(x=age, y=2), shape = 0 , size=2, color="gray80")+
  geom_point(data= data_scheme_02_b, aes(x=age, y= 2), shape = 15, size= 2, color="gray20")+
  geom_segment(data= data_scheme_02_b, aes(x = age, y = 2.1, xend = age, yend = 2.9), color = "gray80", lty=3)+
  geom_point(data= data_scheme_02_b, aes(x=age, y= 3), shape = 15, size= 2, color="green")

dia03_B1

dia03_B2<- ggplot()+
  theme_classic()+
  theme (  axis.text.y = element_blank(),
           axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.text.x = element_blank(),
           axis.title.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.line.x = element_blank(),
           axis.line.y = element_blank())+
  labs(x="Age (cal yr BP)")+
  coord_cartesian(xlim = c(0,8000))+
  geom_point(data= data_scheme_02_b, aes(x=age, y= 1), shape = 15, size= 2, color="green")+
  geom_segment(aes(x = data_scheme_02_b$age[-length(data_scheme_02_b$age)]+100, y = 1.1, xend = data_scheme_02_b$age[-1]-100, yend = 1.1), color = "gray50", size = 3)+
  geom_segment(aes(x = data_scheme_03_b$age, y = 1.4,
                   xend = data_scheme_03_b$age, yend = data_scheme_03_b$VALUE+0.9), color = "gray50", lty = 3)+
  # geom_line(aes(x=data_scheme_03_b$age, y=data_scheme_03_b$VALUE+1), lty=2)+
  geom_point(aes(x=data_scheme_03_b$age, y=data_scheme_03_b$VALUE+1), shape=18, color="green", size=3)

dia03_B2

dia03_C1 <- ggplot()+
  theme_classic()+
  theme (axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.title.x= element_blank())+
  labs(x="Age (cal yr BP)")+
  coord_cartesian(xlim = c(0,8000))+
  geom_point(data = data_scheme_01, aes(x=age, y=1.3), shape = 0 , size=2, color="gray20")+
  geom_segment(aes(x = 4e3, y = 1.4, xend = 4e3, yend = 1.7), arrow = arrow(length = unit(0.3, "cm")), color = "gray80")+
  geom_segment(aes(x = BINs_b, y = 1.75, xend = BINs_b, yend = 2.25), color = "red")+
  geom_point(data = data_scheme_01, aes(x=age, y=2), shape = 0 , size=2, color="gray80")+
  geom_point(data= data_scheme_02_c, aes(x=age, y= 2), shape = 15, size= 2, color="gray20")+
  geom_segment(data= data_scheme_02_c, aes(x = age, y = 2.1, xend = age, yend = 2.9), lty = 3, color = "gray80")+
  geom_point(data= data_scheme_02_c, aes(x=age, y= 3), shape = 15, size= 2, color="red")

dia03_C1

dia03_C2<-ggplot()+
  theme_classic()+
  theme (  axis.text.y = element_blank(),
           axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.text.x = element_blank(),
           axis.title.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.line.x = element_blank(),
           axis.line.y = element_blank())+
  labs(x="Age (cal yr BP)")+
  coord_cartesian(xlim = c(0,8000))+
  geom_point(data= data_scheme_02_c, aes(x=age, y= 1), shape = 15, size= 2, color="red")+
  geom_segment(aes(x = data_scheme_02_c$age[-length(data_scheme_02_c$age)]+100, y = 1.1, xend = data_scheme_02_c$age[-1]-100, yend = 1.1), color = "gray50", size = 3)+
  geom_segment(aes(x = data_scheme_03_c$age, y = 1.4,
                   xend = data_scheme_03_c$age, yend = data_scheme_03_c$VALUE+0.9), 
               lty = 3, color = "gray50")+
  #geom_line(aes(x=data_scheme_03_c$age, y=data_scheme_03_c$VALUE+1), lty=2)+
  geom_point(aes(x=data_scheme_03_c$age, y=data_scheme_03_c$VALUE+1), shape=18, color="red", size=3)


dia03_C2

dia03_sum <- ggarrange(
  dia03_A2, dia03_B2, dia03_C2,
  dia03_A1, dia03_B1, dia03_C1,
  nrow = 2, ncol = 3
)

dia03_sum_anot <- annotate_figure(dia03_sum, bottom = "Age (yr BP)")

dia03_sum_anot

dia03_top <-ggplot()+
  theme_classic()+
  theme (axis.text.y = element_blank(),
         axis.ticks.y = element_blank())+
  labs(x="Age (cal yr BP)",
       y= "Rate-of-Change score",
       title = "Summarisation results from all Mowing Windows")+
  coord_cartesian(xlim = c(0,8000))+
  geom_point(aes(x=data_scheme_03_a$age, y=data_scheme_03_a$VALUE), shape=18, color="blue", size=3)+
  geom_point(aes(x=data_scheme_03_b$age, y=data_scheme_03_b$VALUE), shape=18, color="green", size=3)+
  geom_point(aes(x=data_scheme_03_c$age, y=data_scheme_03_c$VALUE), shape=18, color="red", size=3)+
  geom_line(data=data_scheme_03, aes(x=age, y=VALUE), color="gray30", lty=2)

dia03 <- ggarrange(
  dia03_top,
  dia03_sum_anot,
  nrow=2, heights = c(0.8,1)
) 

dia_fin <- ggarrange(
  ggarrange(dia01,
            dia02, nrow = 2, heights = c(1,0.5), labels = c("A","B")),
  dia03,
  nrow = 1, labels = c("","C")
)

dia_fin

ggsave("~/RESULTS/Methods/FIN/Supplementary_F6A.pdf",
       plot = dia_fin,
       height = 20, width = 35, units="cm")


data_scheme_03

data_scheme_04 <- matrix (nrow=nrow(data_scheme_03), ncol= 10)

data_scheme_04[,1] <- data_scheme_03$VALUE


for(i in 2:10){
  
  for(k in 1:nrow(data_scheme_04)){
    data_scheme_04[k,i] <- runif(1, 
                                 min = data_scheme_04[k,1]- sd(data_scheme_04[,1]),
                                 max = data_scheme_04[k,1]+ sd(data_scheme_04[,1]))  
  }
  
  data_scheme_04[,i] <- lowess(data_scheme_03$age,data_scheme_04[,i],f=.2,iter=100)$y
}



data_scheme_04 <-as.data.frame(data_scheme_04)

data_scheme_04$AGE <- data_scheme_03$age


dia04 <- ggarrange(
  data_scheme_04 %>%
    pivot_longer(cols = -c(AGE)) %>%
    ggplot(aes(y=value, x= AGE))+
    theme_classic()+
    geom_line(aes(group= name), lty= 2, color="gray50")+
    labs(x="Age (cal yr BP)",
         y= "Rate-of-Change score",
         title ="Summarisation results from all randomisations")+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  ,
  data_scheme_04 %>%
    pivot_longer(cols = -c(AGE)) %>%
    group_by(AGE) %>%
    summarise(VALUE = median (value),
              VALUE.dw = quantile(value, 0.2),
              VALUE.up = quantile(value, 0.8)
              
    ) %>%
    ggplot(aes(y=VALUE, x= AGE))+
    theme_classic()+
    geom_ribbon(aes(ymin =VALUE.dw, ymax=VALUE.up),color="gray80", fill="gray80")+
    geom_line(color="gray30", size = 1)+
    labs(x="Age (cal yr BP)",
         y= "Rate-of-Change score")+
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank())
  ,
  data_scheme_04 %>%
    pivot_longer(cols = -c(AGE)) %>%
    group_by(AGE) %>%
    summarise(VALUE = median (value),
              VALUE.dw = quantile(value, 0.2),
              VALUE.up = quantile(value, 0.8)
              
    ) %>%
    mutate(
      pred.gam = predict.gam(gam(VALUE~s(AGE,k=3), data = .)),
      pred.gam.diff = VALUE - pred.gam,
      Peak = (pred.gam.diff) > 1.5*sd(pred.gam.diff),
      UP = VALUE > pred.gam
    ) %>%
    ggplot(aes(y=VALUE, x= AGE))+
    theme_classic()+
    geom_ribbon(aes(ymin =VALUE.dw, ymax=VALUE.up),color="gray80", fill="gray80")+
    geom_line(color="gray30", size = 1)+
    geom_line(aes(y=pred.gam), color= "blue", size = 2)+
    geom_point(color="gray30")+
    geom_segment(data= . %>% filter(UP == T & Peak == F), aes(y=pred.gam,
                                                              yend=VALUE,
                                                              x=AGE, 
                                                              xend=AGE), lty= 3)+
    geom_segment(data= . %>% filter(UP == T & Peak == T),
                 aes(y=pred.gam,
                     yend=VALUE,
                     x=AGE, 
                     xend=AGE), color= "red", lty = 3, size =1)+
    geom_point(data = . %>%filter(Peak == T), color = "black", size= 6)+
    geom_point(data = . %>%filter(Peak == T), color = "green", size= 5)+
    labs(x="Age (cal yr BP)",
         y= "Rate-of-Change score",
         title = "Detection of Peak points")+
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank())
,  
  nrow = 1
)
dia04

dia04_a <- annotate_figure(dia04,
                            left= "Rate-of-Change score",
                            bottom = "Age (cal yr BP)")
dia04_a

ggsave("~/RESULTS/Methods/FIN/Supplementary_F6B.pdf",
       plot = dia04_a,
       height = 12, width = 20, units="cm")


# bonus 

age.df <- data.frame(AGE = data_site_A$list_ages$ages$age,
       SAMPLE = data_site_A$list_ages$age_position %>%
         as_tibble() %>%
         sample_n(.,1) %>%
         t() %>%
         as_tibble())


data_site_A$list_ages$age_position %>%
  t() %>%
  as_tibble() %>%
  mutate(AGE = data_site_A$list_ages$ages$age) %>%
  pivot_longer(cols = -c(AGE)) %>%
  sample_n(.,10e3) %>% 
  arrange(AGE) %>% 
  ggplot(aes(x=AGE, y=value-AGE))+
  geom_point(alpha= 1, color="gray80", shape = 0)+
  geom_line(data = age.df , aes(x=AGE,y=V1-AGE), color="gray30")+
  geom_point(data = age.df , aes(x=AGE,y=V1-AGE), shape = 15, color = "gray30")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("~/RESULTS/Methods/FIN/Supplementary_F6C.pdf",
       height = 12, width = 20, units="cm")


ggarrange(
  data_site_A$filtered.counts %>%
    as_tibble() %>%
    mutate(SUM = rowSums(.),
           AGE = data_site_A$list_ages$ages$age) %>%
    pivot_longer(cols= -c(AGE, SUM)) %>%
    filter(AGE < 8000) %>%
    mutate(VALUE = value/SUM *100) %>%
    filter(VALUE > 10) %>%
    ggplot(aes(x=AGE, y= value))+
    geom_bar(aes(fill=name), stat = "identity", orientation = "x")+
    geom_hline(yintercept = 150)+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  ,
  data_site_A$filtered.counts %>%
    as_tibble() %>%
    mutate(SUM = rowSums(.),
           AGE = data_site_A$list_ages$ages$age) %>%
    pivot_longer(cols= -c(AGE, SUM)) %>%
    filter(AGE < 8000) %>%
    mutate(VALUE = value/SUM *100) %>%
    filter(VALUE > 10) %>%
    ggplot(aes(x=AGE, y= value))+
    geom_bar(aes(fill=name), stat = "identity", orientation = "x", position= position_fill())+
    theme_classic()+
        theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
)


ggsave("~/RESULTS/Methods/FIN/Supplementary_F6D.pdf",
       height = 7, width = 20, units="cm")


####################################
#               SAVE               #
####################################

# save.image("C:/Users/omo084/OneDrive - University of Bergen/PRIVATE/ROC_method/ENV_METHOD_2020807.RData")

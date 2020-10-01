# load("C:/Users/omo084/OneDrive - University of Bergen/PRIVATE/ROC_method/ENV_METHOD_2020807.RData")

# ----------------------------------------------
#                     SETUP
# ----------------------------------------------
library(tidyverse)
library(maps)
library(RColorBrewer)
library(MuMIn)
library(glmmTMB)
library(emmeans)
library(parallel)
library(doParallel)
library (RColorBrewer)
library(cowplot)
library(ggpubr)
library(viridis)
library(devtools)
devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package")

library(RRatepol)

theme_set(theme_classic())

# ----------------------------------------------
#           LOAD DATA AND FUNCTIONS
# ----------------------------------------------
data <- RRatepol::example_data

files.sources <- list.files("functions/") 
sapply(paste0("functions/", files.sources, sep =""), source)

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
time_seq <- data$list_ages[[4]]$ages$age

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

perform_sim_ld_recent_levels <- fc_test_simlutated_data_succsess(sim_ld_recent_levels, breaks = breaks_recent)
perform_sim_ld_late_levels <- fc_test_simlutated_data_succsess(sim_ld_late_levels, breaks = breaks_late)
perform_sim_hd_recent_levels <- fc_test_simlutated_data_succsess(sim_hd_recent_levels, breaks = breaks_recent)
perform_sim_hd_late_levels <- fc_test_simlutated_data_succsess(sim_hd_late_levels, breaks = breaks_late)

perform_sim_ld_recent_BINs <- fc_test_simlutated_data_succsess(sim_ld_recent_BINs, breaks = breaks_recent)
perform_sim_ld_late_BINs <- fc_test_simlutated_data_succsess(sim_ld_late_BINs, breaks = breaks_late)
perform_sim_hd_recent_BINs <- fc_test_simlutated_data_succsess(sim_hd_recent_BINs, breaks = breaks_recent)
perform_sim_hd_late_BINs <- fc_test_simlutated_data_succsess(sim_hd_late_BINs, breaks = breaks_late)

perform_sim_ld_recent_MW <- fc_test_simlutated_data_succsess(sim_ld_recent_MW, breaks = breaks_recent)
perform_sim_ld_late_MW <- fc_test_simlutated_data_succsess(sim_ld_late_MW, breaks = breaks_late)
perform_sim_hd_recent_MW <- fc_test_simlutated_data_succsess(sim_hd_recent_MW, breaks = breaks_recent)
perform_sim_hd_late_MW <- fc_test_simlutated_data_succsess(sim_hd_late_MW, breaks = breaks_late)

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

data_success_sum <- within(data_success_sum, PEAK <- factor(PEAK, levels = c("PEAK.G","PEAK.T","PEAK.S")))

# Produce FIG S2 !!!

# cluster setup

nrCores <- parallel::detectCores()

## FOCUS

data_success_focus <- data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  ungroup() %>%
  dplyr::select(-c(SEGMENT))

mod_success_focus <-  glmmTMB(VALUE~WU*PEAK*Position*Diversity+(WU|dataset.ID),
                           data=data_success_focus,
                           family=betabinomial(link = "probit"))

mod_success_focus <-  glmmTMB(VALUE~PEAK+(WU|dataset.ID),
                              data=data_success_focus,
                              family=betabinomial(link = "probit"))


cl <- parallel::makeCluster(nrCores-1)
doParallel::registerDoParallel(cl); 
parallel::clusterExport(cl,c("data_success_focus","nrCores"),envir=environment());
parallel::clusterEvalQ(cl,library("glmmTMB"))

mod_success_focus_dd <- pdredge(mod_success_focus,cluster = cl, trace = T)

(mod_success_focus_dd$AICc - mod_success_focus_dd$AICc[1])<2

mod_success_focus_dd[1,]
#  -> FULL model 

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
                                   family=betabinomial(link = "probit"))

cl <- parallel::makeCluster(nrCores-1)
doParallel::registerDoParallel(cl); 
parallel::clusterExport(cl,c("data_success_focus_MW_G","nrCores"),envir=environment());
parallel::clusterEvalQ(cl,library("glmmTMB"))

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
                                      family=betabinomial(link = "probit"),
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
                              family=betabinomial(link = "probit"))

cl <- parallel::makeCluster(nrCores-1)
doParallel::registerDoParallel(cl); 
parallel::clusterExport(cl,c("data_success_false","nrCores"),envir=environment());
parallel::clusterEvalQ(cl,library("glmmTMB"))

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
                              family=betabinomial(link = "probit"),
                              control = glmmTMBControl(parallel = nrCores-1))

data_success_false_MW_G <- data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  filter(WU == "MW") %>% 
  filter(PEAK == "PEAK.G") %>%
  ungroup() %>%
  dplyr::select(-c(SEGMENT,WU,PEAK))

mod_success_false_MW_G <-  glmmTMB(VALUE~Position*Diversity*DC*SMOOTH+(1|dataset.ID),
                                   data=data_success_false_MW_G,
                                   family=betabinomial(link = "probit"))


cl <- parallel::makeCluster(nrCores-1)
doParallel::registerDoParallel(cl); 
parallel::clusterExport(cl,c("data_success_false_MW_G","nrCores"),envir=environment());
parallel::clusterEvalQ(cl,library("glmmTMB"))

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
                                      family=betabinomial(link = "probit"),
                                      control = glmmTMBControl(parallel = nrCores-1))


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

data_site_A <- list(dataset.id = data$dataset.id[[4]],
                    filtered.counts = data$filtered.counts[[4]],
                    list_ages = data$list_ages[[4]])

data_site_A_dom <- get_dominant_pollen_taxa(data_site_A)

data_site_A$filtered.counts %>%
  as_tibble() %>%
  dim()

#Site B

data_site_B <- list(dataset.id = data$dataset.id[[1]],
                    filtered.counts = data$filtered.counts[[1]],
                    list_ages = data$list_ages[[1]])

data_site_B_dom <- get_dominant_pollen_taxa(data_site_B)

data_site_B$filtered.counts %>%
  as_tibble() %>%
  dim()


# SITE c

data_site_C <- list(dataset.id = data$dataset.id[[2]],
                    filtered.counts = data$filtered.counts[[2]],
                    list_ages = data$list_ages[[2]])

data_site_C_dom <- get_dominant_pollen_taxa(data_site_C)

data_site_C$filtered.counts %>%
  as_tibble() %>%
  dim()


# Site D

data_site_D <- list(dataset.id = data$dataset.id[[3]],
                    filtered.counts = data$filtered.counts[[3]],
                    list_ages = data$list_ages[[3]])

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


getPalette = colorRampPalette(brewer.pal(9, "Set1"))
Palette.1<- getPalette(length(common_taxa))
names(Palette.1)<- sort(common_taxa)  


fc_get_pollen_data <- function (data, smooth_method , Common.list)
{
  # remove the sample ID
  if (is.numeric(unlist(data$filtered.counts[,1]))==F){
    data$filtered.counts <- data$filtered.counts[,-1]
  }
  
  data.ext <- RRatepol::fc_extract_data(data$filtered.counts,
                                        data$list_ages) %>%
    RRatepol::fc_smooth_pollen_data(.,smooth_method  = smooth_method ,
                                    smooth_N_points  = 5,
                                    smooth_N_max  = 9,
                                    smooth_age_range  = 500) %>%
    RRatepol::fc_check_data(.,proportion = T)
  
  plot.data <- data.ext$Pollen %>%
    select(any_of(Common.list)) %>%
    rownames_to_column() %>%
    pivot_longer(cols = -c(rowname)) %>%
    rename(sample.id = rowname) %>%
    inner_join(.,data.ext$Age, by="sample.id")
  
  return (plot.data)
}


data_site_A_pollen <-fc_get_pollen_data(data_site_A, smooth_method  = "none",common_taxa)

data_site_B_pollen <-fc_get_pollen_data(data_site_B, smooth_method  = "none", common_taxa)

data_site_C_pollen <-fc_get_pollen_data(data_site_C, smooth_method  = "none",common_taxa)

data_site_D_pollen <-fc_get_pollen_data(data_site_D, smooth_method  = "none",common_taxa)


# Rate of Change

data_site_A_RoC_levels <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_A$filtered.counts,
                                                    data_source_age = data_site_A$list_ages,
                                                    smooth_method  = "age.w", 
                                                    smooth_N_points  = 5,
                                                    smooth_age_range  = 500, 
                                                    smooth_N_max   = 9,
                                                    Working_Units  = "levels",
                                                    bin_size  = 500,
                                                    Number_of_shifts = 5,
                                                    rand = 10000,
                                                    standardise = T, 
                                                    N_pollen_grains  = 150, 
                                                    DC = "chisq",
                                                    interest_threshold  = age_lim,
                                                    Debug = F)

data_site_A_RoC_BINs <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_A$filtered.counts,
                                                  data_source_age = data_site_A$list_ages,
                                                  smooth_method  = "age.w", 
                                                  smooth_N_points  = 5,
                                                  smooth_age_range = 500, 
                                                  smooth_N_max   = 9,
                                                  Working_Units  = "BINs",
                                                  bin_size  = 500,
                                                  Number_of_shifts  = 5,
                                                  rand = 10000,
                                                  standardise = T, 
                                                  N_pollen_grains  = 150, 
                                                  DC = "chisq",
                                                  interest_threshold  = age_lim,
                                                  Debug = F)

data_site_A_RoC_MW <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_A$filtered.counts,
                                                data_source_age = data_site_A$list_ages,
                                                smooth_method  = "age.w", 
                                                smooth_N_points  = 5,
                                                smooth_age_range = 500, 
                                                smooth_N_max   = 9,
                                                Working_Units  = "MW",
                                                bin_size  = 500,
                                                Number_of_shifts  = 5,
                                                rand = 10000,
                                                standardise = T, 
                                                N_pollen_grains  = 150, 
                                                DC = "chisq",
                                                interest_threshold  = age_lim,
                                                Debug = F)

data_Site_B_RoC_levels <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_B$filtered.counts,
                                                    data_source_age = data_site_B$list_ages,
                                                    smooth_method  = "shep", 
                                                    smooth_N_points  = 5,
                                                    smooth_age_range = 500, 
                                                    smooth_N_max   = 9,
                                                    Working_Units  = "levels",
                                                    bin_size  = 500,
                                                    Number_of_shifts  = 5,
                                                    rand = 10000,
                                                    standardise = T, 
                                                    N_pollen_grains  = 150, 
                                                    DC = "chord",
                                                    interest_threshold  = age_lim,
                                                    Debug = F)

data_Site_B_RoC_BINs <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_B$filtered.counts,
                                                  data_source_age = data_site_B$list_ages,
                                                  smooth_method  = "shep", 
                                                  smooth_N_points  = 5,
                                                  smooth_age_range = 500, 
                                                  smooth_N_max   = 9,
                                                  Working_Units  = "BINs",
                                                  bin_size  = 500,
                                                  Number_of_shifts  = 5,
                                                  rand = 10000,
                                                  standardise = T, 
                                                  N_pollen_grains  = 150, 
                                                  DC = "chord",
                                                  interest_threshold  = age_lim,
                                                  Debug = F)

data_Site_B_RoC_MW <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_B$filtered.counts,
                                                data_source_age = data_site_B$list_ages,
                                                smooth_method  = "shep", 
                                                smooth_N_points  = 5,
                                                smooth_age_range = 500, 
                                                smooth_N_max   = 9,
                                                Working_Units  = "MW",
                                                bin_size  = 500,
                                                Number_of_shifts  = 5,
                                                rand = 10000,
                                                standardise = T, 
                                                N_pollen_grains  = 150, 
                                                DC = "chord",
                                                interest_threshold  = age_lim,
                                                Debug = F)

data_Site_C_RoC_levels <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_C$filtered.counts,
                                                    data_source_age = data_site_C$list_ages,
                                                    smooth_method  = "shep", 
                                                    smooth_N_points  = 5,
                                                    smooth_age_range = 500, 
                                                    smooth_N_max   = 9,
                                                    Working_Units  = "levels",
                                                    bin_size  = 500,
                                                    Number_of_shifts  = 5,
                                                    rand = 10000,
                                                    standardise = T, 
                                                    N_pollen_grains  = 150, 
                                                    DC = "chord",
                                                    interest_threshold  = age_lim,
                                                    Debug = F)

data_Site_C_RoC_BINs <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_C$filtered.counts,
                                                  data_source_age = data_site_C$list_ages,
                                                  smooth_method  = "shep", 
                                                  smooth_N_points  = 5,
                                                  smooth_age_range = 500, 
                                                  smooth_N_max   = 9,
                                                  Working_Units  = "BINs",
                                                  bin_size  = 500,
                                                  Number_of_shifts  = 5,
                                                  rand = 10000,
                                                  standardise = T, 
                                                  N_pollen_grains  = 150, 
                                                  DC = "chord",
                                                  interest_threshold  = age_lim,
                                                  Debug = F)

data_Site_C_RoC_MW <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_C$filtered.counts,
                                                data_source_age = data_site_C$list_ages,
                                                smooth_method  = "shep", 
                                                smooth_N_points  = 5,
                                                smooth_age_range = 500, 
                                                smooth_N_max   = 9,
                                                Working_Units  = "MW",
                                                bin_size  = 500,
                                                Number_of_shifts  = 5,
                                                rand = 10000,
                                                standardise = T, 
                                                N_pollen_grains  = 150, 
                                                DC = "chord",
                                                interest_threshold  = age_lim,
                                                Debug = F)

data_Site_D_RoC_levels <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_D$filtered.counts,
                                                    data_source_age = data_site_D$list_ages,
                                                    smooth_method  = "shep", 
                                                    smooth_N_points  = 5,
                                                    smooth_age_range = 500, 
                                                    smooth_N_max   = 9,
                                                    Working_Units  = "levels",
                                                    bin_size  = 500,
                                                    Number_of_shifts  = 5,
                                                    rand = 10000,
                                                    standardise = T, 
                                                    N_pollen_grains  = 150, 
                                                    DC = "chord",
                                                    interest_threshold  = age_lim,
                                                    Debug = F)

data_Site_D_RoC_BINs <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_D$filtered.counts,
                                                  data_source_age = data_site_D$list_ages,
                                                  smooth_method  = "shep", 
                                                  smooth_N_points  = 5,
                                                  smooth_age_range = 500, 
                                                  smooth_N_max   = 9,
                                                  Working_Units  = "BINs",
                                                  bin_size  = 500,
                                                  Number_of_shifts  = 5,
                                                  rand = 10000,
                                                  standardise = T, 
                                                  N_pollen_grains  = 150, 
                                                  DC = "chord",
                                                  interest_threshold  = age_lim,
                                                  Debug = F)


data_Site_D_RoC_MW <- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_D$filtered.counts,
                                                data_source_age = data_site_D$list_ages,
                                                smooth_method  = "shep", 
                                                smooth_N_points  = 5,
                                                smooth_age_range = 500, 
                                                smooth_N_max   = 9,
                                                Working_Units  = "MW",
                                                bin_size  = 500,
                                                Number_of_shifts  = 5,
                                                rand = 10000,
                                                standardise = T, 
                                                N_pollen_grains  = 150, 
                                                DC = "chord",
                                                interest_threshold  = age_lim,
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
  )+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))



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
    y="Proportion of total pollen"
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
    y="Proportion of total pollen"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        strip.background = element_blank(),
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
    y="Proportion of total pollen"
  )+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        strip.background = element_blank(),
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
    y="Proportion of total pollen"
  )+
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y  = element_blank(),
        strip.background = element_blank(),
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
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y  = element_blank(),
  )

FIG3_Site_D_3_BINs <-  data_Site_D_RoC_BINs %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
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
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y  = element_blank(),
  )
FIG3_Site_D_3_MW <-  data_Site_D_RoC_MW %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,2))+
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
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
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
  labels = c("A","B","C","D"),
  heights = c(1.8,1,1,1.2),
  ncol = 1, nrow = 4, legend = "none")
FIG3_Site_comparison



ggsave("METHOD_RESULTS//FIG3_v02.pdf",
       width = 20, height = 22, units = "cm")





############
#   FIG 4
############

fc_calculate_RoC_comparison <- function(data, Working_Units , bin_size , Number_of_shifts , rand ,peek, interest_threshold ){
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  
  for(i in 1:20)
  {
    data.temp<- RRatepol::fc_estimate_RoC( data_source_pollen =  data$filtered.counts,
                              data_source_age = data$list_ages,
                              smooth_method  = performance.smooth[i],
                              smooth_N_points  = 5,
                              smooth_age_range = 500, 
                              smooth_N_max   = 9,
                              Working_Units  = Working_Units ,
                              bin_size  = bin_size ,
                              Number_of_shifts  = Number_of_shifts ,
                              rand = rand,
                              standardise = F, 
                              N_pollen_grains  = 150 ,
                              DC = performance.DC[i],
                              interest_threshold  = interest_threshold ,
                              Peak = peek,
                              Debug = F) %>%
      as_tibble()
    
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


data_example_MW <- fc_calculate_RoC_comparison(data_site_A,
                                               Working_Units  = "MW",
                                               bin_size  = 500,
                                               Number_of_shifts  = 5,
                                               rand = 10000,
                                               peek = "GAM",
                                               interest_threshold  =  age_lim)


data_example_MW <- within(data_example_MW, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
levels(data_example_MW$DC) <- c("Euc","Euc.sd","Chord","Chisq")

data_example_MW <- within(data_example_MW, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))
levels(data_example_MW$SMOOTH)<-c("None","M.avg","Grimm","Age.w","Shep")

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


###############
# SUPLEMENTARY 
###############

##########
# FIG S1 #
########## 

# Example of simulation of enviroemntal data

genererate_exmaple_data <- function(time= time_seq,
                                    nforc=4,
                                    mean=100,
                                    sdev=.15, 
                                    nprox=50, 
                                    var=20,
                                    range=15,
                                    manual.edit = T,
                                    breaks,
                                    jitter = T,
                                    rarity=T) {
  
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
  return(list(forcing = forcing, pollen = data.source.pollen, age = data.source.age))
}

Supplementary_F1a_data <- genererate_exmaple_data(breaks = c(2000,3000))

Supplementary_F1a <-ggarrange(as.data.frame(Supplementary_F1a_data$forcing) %>%
                               mutate(AGE = time) %>%
                               pivot_longer(., cols = -c(AGE)) %>%
                               arrange(AGE,value) %>%
                               ggplot(aes(x=AGE, y= value))+
                               geom_vline(xintercept = c(2000,3000), color="gray80", size=0.1)+
                               geom_line(aes(color=name))+
                               theme_classic()+
                               coord_flip(xlim=c(8000,0))+
                               scale_x_continuous(trans = "reverse")+
                               theme(
                                 axis.ticks.y = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.title.y = element_blank(),
                                     legend.position = "none"
                                     )+
                                ylab("Value of env. variable")+
                                xlab("Age (cal yr BP)"),
                             fc_extract_data(Supplementary_F1a_data$pollen, Supplementary_F1a_data$age) %>%
                               fc_smooth_pollen_data("none") %>%
                               fc_check_data(., proportion = T) %>%
                               pluck("Pollen") %>%
                               mutate(AGE = time) %>%
                               pivot_longer(-c(AGE)) %>%
                               ggplot(aes(x=AGE, y=value))+
                               geom_ribbon(aes(ymin=rep(0,length(value)), ymax=value, fill=name), 
                                           color="gray20", alpha=1/5, size=0.1)+
                               geom_vline(xintercept = c(2000,3000), color="gray80", size=0.1)+
                               coord_flip(xlim=c(8000,0), ylim=c(0,1))+
                               theme_classic()+
                               scale_x_continuous(trans = "reverse")+
                               xlab("")+
                               ylab("Proportion of pollen grains")+
                               theme(
                                     axis.ticks.y = element_blank(),
                                     axis.text.y = element_blank(),
                                     legend.position = "none"),
                             ncol=2, align = "h")


Supplementary_F1a

Supplementary_F1b_data <- genererate_exmaple_data(breaks = c(5500,6500))

Supplementary_F1b <-ggarrange(as.data.frame(Supplementary_F1b_data$forcing) %>%
                                mutate(AGE = time) %>%
                                pivot_longer(., cols = -c(AGE)) %>%
                                arrange(AGE,value) %>%
                                ggplot(aes(x=AGE, y= value))+
                                geom_vline(xintercept = c(5500,6500), color="gray80", size=0.1)+
                                geom_line(aes(color=name))+
                                theme_classic()+
                                coord_flip(xlim=c(8000,0))+
                                scale_x_continuous(trans = "reverse")+
                                theme(
                                  axis.ticks.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.title.y = element_blank(),
                                  legend.position = "none"
                                )+
                                ylab("Value of env. variable")+
                                xlab("Age (cal yr BP)"),
                              fc_extract_data(Supplementary_F1b_data$pollen, Supplementary_F1b_data$age) %>%
                                fc_smooth_pollen_data("none") %>%
                                fc_check_data(., proportion = T) %>%
                                pluck("Pollen") %>%
                                mutate(AGE = time) %>%
                                pivot_longer(-c(AGE)) %>%
                                ggplot(aes(x=AGE, y=value))+
                                geom_ribbon(aes(ymin=rep(0,length(value)), ymax=value, fill=name), 
                                            color="gray20", alpha=1/5, size=0.1)+
                                geom_vline(xintercept = c(5500,6500), color="gray80", size=0.1)+
                                coord_flip(xlim=c(8000,0), ylim=c(0,1))+
                                theme_classic()+
                                scale_x_continuous(trans = "reverse")+
                                xlab("")+
                                ylab("Proportion of pollen grains")+
                                theme(
                                  axis.ticks.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  legend.position = "none"),
                              ncol=2, align = "h")


Supplementary_F1b

Supplementary_F1 <- ggarrange(
  Supplementary_F1a,Supplementary_F1b,
  nrow = 2, labels=c("A","B")
)


Supplementary_F1_fin <- annotate_figure(Supplementary_F1, left= "Age (cal yr BP)")

Supplementary_F1_fin

ggsave("METHOD_RESULTS/Supplementary_F1.pdf",
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
            coord_cartesian(ylim = c(0,0.4))+
            scale_color_manual(values = Color.legen_Diversity )+
            scale_fill_manual(values = Color.legen_Diversity )+
            labs(y="",x=""),
          nrow = 1, common.legend = T, legend = "none") ,
nrow = 4) %>% annotate_figure(left = "Proportion of peak detection")


ggsave("METHOD_RESULTS//Supplementary_F2.pdf",
       width = 20, height = 20, units = "cm")





##########
# FIG S3 #
########## 

for(i in seq(100,1e3,100)){
 df_w<- RRatepol::fc_estimate_RoC(data_source_pollen = data_site_A$filtered.counts,
                            data_source_age = data_site_A$list_ages,
                            smooth_method  = "age.w", 
                            smooth_N_points  = 5,
                            smooth_age_range = 500, 
                            smooth_N_max   = 9,
                            Working_Units  = "MW",
                            bin_size  = i,
                            Number_of_shifts  = 5,
                            rand = 1000,
                            standardise = T, 
                            N_pollen_grains  = 150, 
                            DC = "chisq",
                            interest_threshold  = age_lim,
                            Debug = F)
    if(exists("BIM_tibble")==F){
      BIM_tibble <- tibble(df_w, BIN =i)
    } else {
      BIM_tibble <- bind_rows(BIM_tibble,tibble(df_w, BIN =i))
    }
}


Supplementary_F5 <- BIM_tibble %>%
  ggplot()+
  geom_rug(data = data_site_A$list_ages$ages, aes(x=age), sides = "b")+
  geom_line(aes(x= AGE, y= ROC, group=BIN, color= BIN))+
  viridis::scale_colour_viridis(direction = -1)+
  labs(x="Age (cal yr BP)",
       y= "Rate-of-Change score",
       color = "Time bin size")



ggsave("METHOD_RESULTS/Supplementary_F3.pdf",
       plot = Supplementary_F5,
       height = 12, width = 20, units="cm")





####################################
#               SAVE               #
####################################

# save.image("C:/Users/omo084/OneDrive - University of Bergen/PRIVATE/ROC_method/ENV_METHOD_2020807.RData")

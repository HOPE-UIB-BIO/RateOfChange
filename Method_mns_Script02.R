# load("~/DATA/temp/ENV_METHOD_20200525.RData")

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

# ----------------------------------------------
#             LOAD DATA & FUNCTIONS
# ----------------------------------------------
# download.file("https://www.dropbox.com/s/3hp7rv03mkg4pjz/tibble_Europe_filtered05.03.20.RData?dl=1","~/input/DATA/tibble_Europe_filtered05.03.20.RData")

setwd("~/GITHUB/RateOfChange")

load("~/DATA/input/tibble_Europe_filtered05.03.20.RData")

files.sources <- list.files("~/GITHUB/RateOfChange/functions/") 
sapply(paste0("~/GITHUB/RateOfChange/functions/", files.sources, sep =""), source)

glimpse(tibble_Europe2)

fc_calculate_RoC_comparison <- function(data, Working.Unit, BIN.size, N.shifts, rand ,interest.treshold){
  
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
                            Debug = F) %>%
      as_tibble()
    
    # PEAK detection 
    # Median peak treshold
    # treshold for RoC peaks is set as median of all RoC in dataset
    r.treshold <- median(data.temp$ROC)
    # mark peaks which have 95% quantile above the treshold asPeak.treshold
    data.temp$PEAK.T <- data.temp$ROC.dw > r.treshold
    
    # GAM  
    # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
    pred.gam <-  predict.gam(gam(ROC~s(AGE), data = data.temp))
    pred.gam.diff <- data.temp$ROC - pred.gam
    data.temp$PEAK.G <- (pred.gam.diff) > 1.5*sd(pred.gam.diff)
    
    # SNI  
    # set moving window of 5 times higher than average distance between samples
    mean.age.window <- 5 * mean( diff(data.temp$AGE) )
    # calculate SNI (singal to noise ratio)
    SNI.calc <- CharSNI(data.frame(data.temp$AGE, data.temp$ROC, pred.gam),mean.age.window)
    # mark points with SNI higher than 3
    data.temp$PEAK.S <- SNI.calc$SNI > 3 & data.temp$ROC > pred.gam
    
    data.temp.sum<-data.temp %>%
      mutate(SMOOTH = performance.smooth[i],
             DC = performance.DC[i]) %>%
      select(SMOOTH, DC, Working_Unit ,AGE, ROC,ROC.up, ROC.dw, PEAK.T, PEAK.G, PEAK.S)
    
    if(i == 1){
      res.tibble <- data.temp.sum} else {
        res.tibble <- rbind(res.tibble,data.temp.sum)
      }
  }
  return(res.tibble)
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
time_seq <- tibble_Europe2$list_ages[[2]]$ages$age

# maximal time of focus 
age_lim <- 8000


# -----------------------------------------
#
#             STATISTICAL COMP
# 
# -----------------------------------------

# data simulation

sim_ld_recent_MW <- fc_simulate_pollen_data_in_all_methods(time=time_seq, 
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
                                                           Working.Unit="MW", 
                                                           BIN.size=500, 
                                                           N.shifts=5, 
                                                           N.datasets=N_rep, 
                                                           interest.treshold=8000)

sim_ld_late_MW <- fc_simulate_pollen_data_in_all_methods(time=time_seq, 
                                                           nforc=N_env, 
                                                           mean=100, 
                                                           sdev=.15, 
                                                           nprox=low_diversity, 
                                                           var=20, 
                                                           range=15,
                                                           manual.edit = T,
                                                           breaks=breaks_late,
                                                           jitter = T,
                                                           rarity=T,
                                                           Working.Unit="MW", 
                                                           BIN.size=500, 
                                                           N.shifts=5, 
                                                           N.datasets=N_rep, 
                                                           interest.treshold=8000)

sim_hd_recent_MW <- fc_simulate_pollen_data_in_all_methods(time=time_seq, 
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
                                                           Working.Unit="MW", 
                                                           BIN.size=500, 
                                                           N.shifts=5, 
                                                           N.datasets=N_rep, 
                                                           interest.treshold=8000)

sim_hd_late_MW <- fc_simulate_pollen_data_in_all_methods(time=time_seq, 
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
                                                         Working.Unit="MW", 
                                                         BIN.size=500, 
                                                         N.shifts=5, 
                                                         N.datasets=N_rep, 
                                                         interest.treshold=8000)

# -----------------------------------------
#
#       SUCCESS COMPARISON
# 
# -----------------------------------------

perform_sim_ld_recent_MW<- fc_test_simlutated_data_succsess(sim_ld_recent_MW, breaks = breaks_recent)

perform_sim_ld_late_MW<- fc_test_simlutated_data_succsess(sim_ld_late_MW, breaks = breaks_late)

perform_sim_hd_recent_MW<- fc_test_simlutated_data_succsess(sim_hd_recent_MW, breaks = breaks_recent)

perform_sim_hd_late_MW<- fc_test_simlutated_data_succsess(sim_hd_late_MW, breaks = breaks_late)


data_success_sum <- rbind(
  data.frame(perform_sim_ld_recent_MW$RawData,Position="recent", Diversity= "low"),
  data.frame(perform_sim_ld_late_MW$RawData,Position="late", Diversity= "low"),
  data.frame(perform_sim_hd_recent_MW$RawData,Position="recent", Diversity= "high"),
  data.frame(perform_sim_hd_late_MW$RawData,Position="late", Diversity= "high")
) %>%
  as_tibble()

data_success_sum <- within(data_success_sum, Position <- factor(Position, levels = c("recent","late")))
levels(data_success_sum$Position) <- c("high level density","low level density")

data_success_sum <- within(data_success_sum, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))

data_success_sum <- within(data_success_sum, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))

data_success_sum <- within(data_success_sum, Diversity <- factor(Diversity, levels = c("low","high")))
levels(data_success_sum$Diversity) <- c("low diversity","high diversity")



#data_success_sum$VALUE.S <- data_success_sum$VALUE
#data_success_sum$VALUE.S[data_success_sum$SEGMENT=="empty"] <- 1-data_success_sum$VALUE[data_success_sum$SEGMENT=="empty"]


compare_models <- function(MODEL,SCOPE, order = 1){
  
  res.tib <- tibble(X= SCOPE, VAR =NA,.rows = length(SCOPE) )
  
  for (i in 1: length(SCOPE)){
    
    if(order == 1){
      FORMULA <- paste(". ~ ",SCOPE[i])
    } 
    
    if(order == 2){
      FORMULA <- paste(". ~ . +",SCOPE[i])
    }
    
    new.model <- update(MODEL, FORMULA)  
    new.model.sum <- summary(new.model)
    
    res.tib$VAR[i] <- 1 - new.model.sum$deviance/new.model.sum$null.deviance 
  }
  
  return(res.tib %>% arrange(-VAR))
}


## FOCUS

data_success_sum_G <- data_success_sum %>%
  filter(PEAK == "PEAK.G") %>%
  filter(SEGMENT == "focus") %>%
  ungroup() %>%
  dplyr::select(-c(PEAK,VALUE.S,SEGMENT))


# run 0

scope <- names(data_success_sum_G %>% dplyr::select(-c(dataset.ID,VALUE)))

glm.0 <- glm(VALUE ~ (+1), family = "quasibinomial" , data = data_success_sum_G)

compare_models(glm.0,scope, order = 1)

glm.1 <- glm(VALUE ~ Position, family = "quasibinomial" , data = data_success_sum_G)

anova(glm.0,glm.1, test= "F")


# run 1

scope <- names(data_success_sum_G %>% dplyr::select(-c(dataset.ID,VALUE, Position)))

compare_models(glm.1,scope, order = 2)

glm.2a <- glm(VALUE ~ Position + SMOOTH, family = "quasibinomial" , data = data_success_sum_G)
glm.2b <- glm(VALUE ~ Position * SMOOTH, family = "quasibinomial" , data = data_success_sum_G)

anova(glm.1,glm.2a, test= "F")
anova(glm.2a,glm.2b, test= "F")

# run 2 

scope <- names(data_success_sum_G %>% dplyr::select(-c(dataset.ID,VALUE, Position, SMOOTH)))

compare_models(glm.2b,scope, order = 2)

glm.3a <- glm(VALUE ~ Position * SMOOTH + DC, family = "quasibinomial" , data = data_success_sum_G)
glm.3b <- glm(VALUE ~ Position * SMOOTH * DC, family = "quasibinomial" , data = data_success_sum_G)

anova(glm.2b,glm.3a, test= "F")

# run 3

glm.4a <- glm(VALUE ~ Position * SMOOTH + Diversity , family = "quasibinomial" , data = data_success_sum_G)
glm.4b <- glm(VALUE ~ Position * SMOOTH * Diversity , family = "quasibinomial" , data = data_success_sum_G)

anova(glm.2b,glm.4a, test= "F")
anova(glm.4a,glm.4b, test= "F")

# run 3

glm.5a <- glm(VALUE ~ Position * SMOOTH * Diversity + DC, family = "quasibinomial" , data = data_success_sum_G)
glm.5b <- glm(VALUE ~ Position * SMOOTH * Diversity * DC , family = "quasibinomial" , data = data_success_sum_G)

anova(glm.4b,glm.5a, test= "F")


# GLM FIN

glm.fin <- glm.4b

formula(glm.fin)


## EMPTY

data_success_sum_G <- data_success_sum %>%
  filter(PEAK == "PEAK.G") %>%
  filter(SEGMENT == "empty") %>%
  ungroup() %>%
  dplyr::select(-c(PEAK,VALUE.S,SEGMENT))


# run 0

scope <- names(data_success_sum_G %>% dplyr::select(-c(dataset.ID,VALUE)))

glm.0 <- glm(VALUE ~ (+1), family = "quasibinomial" , data = data_success_sum_G)

compare_models(glm.0,scope, order = 1)

glm.1 <- glm(VALUE ~ Position, family = "quasibinomial" , data = data_success_sum_G)

anova(glm.0,glm.1, test= "F")

# run 1

scope <- names(data_success_sum_G %>% dplyr::select(-c(dataset.ID,VALUE, Position)))

compare_models(glm.1,scope, order = 2)

glm.2a <- glm(VALUE ~ Position + SMOOTH, family = "quasibinomial" , data = data_success_sum_G)
glm.2b <- glm(VALUE ~ Position * SMOOTH, family = "quasibinomial" , data = data_success_sum_G)

anova(glm.1,glm.2a, test= "F")
anova(glm.2a,glm.2b, test= "F")

# run 2 

scope <- names(data_success_sum_G %>% dplyr::select(-c(dataset.ID,VALUE, Position, SMOOTH)))

compare_models(glm.2a,scope, order = 2)

glm.3a <- glm(VALUE ~ Position + SMOOTH + Diversity, family = "quasibinomial" , data = data_success_sum_G)
glm.3b <- glm(VALUE ~ (Position + SMOOTH) * Diversity, family = "quasibinomial" , data = data_success_sum_G)

anova(glm.2a,glm.3a, test= "F")
anova(glm.3a,glm.3b, test= "F")

# run 3

glm.4a <- glm(VALUE ~ (Position + SMOOTH) * Diversity + DC, family = "quasibinomial" , data = data_success_sum_G)
glm.4b <- glm(VALUE ~ (Position + SMOOTH) * Diversity * DC, family = "quasibinomial" , data = data_success_sum_G)

anova(glm.3b,glm.4a, test= "F")
anova(glm.4a,glm.4b, test= "F")

glm.fin.e <- glm.4b

formula(glm.fin.e)

# -----------------------------------------
#
#       MAGNITUDE COMPARISON
# 
# -----------------------------------------

mag_sim_ld_recent_MW <- fc_test_simlutated_data_magnitude(sim_ld_recent_MW)

mag_sim_ld_late_MW <- fc_test_simlutated_data_magnitude(sim_ld_late_MW)

mag_sim_hd_recent_MW <- fc_test_simlutated_data_magnitude(sim_hd_recent_MW)

mag_sim_hd_late_MW <- fc_test_simlutated_data_magnitude(sim_hd_late_MW)


data_mag_sum <- rbind(
  data.frame(mag_sim_ld_recent_MW,Position="recent",Diversity="low"),
  data.frame(mag_sim_ld_late_MW,Position="late",Diversity="low"),
  data.frame(mag_sim_hd_recent_MW,Position="recent",Diversity="high"),
  data.frame(mag_sim_hd_late_MW,Position="late",Diversity="high")
) %>% as_tibble()


# run 0

scope <- names(data_mag_sum %>% dplyr::select(-c(dataset.ID,RoC_max)))

M.glm.0 <- glm(RoC_max ~ (+1), family = "gaussian" , data = data_mag_sum)

compare_models(M.glm.0,scope, order = 1)

M.glm.1 <- glm(RoC_max ~ DC, family = "gaussian" , data = data_mag_sum)

anova(M.glm.0,M.glm.1, test= "Chisq")

# run 1

scope <- names(data_mag_sum %>% dplyr::select(-c(dataset.ID,RoC_max,DC)))

compare_models(M.glm.1,scope, order = 2)

M.glm.2a <- glm(RoC_max ~ DC + Diversity, family = "gaussian" , data = data_mag_sum)
M.glm.2b <- glm(RoC_max ~ DC * Diversity, family = "gaussian" , data = data_mag_sum)

anova(M.glm.1,M.glm.2a, test= "Chisq")
anova(M.glm.2a,M.glm.2b, test= "Chisq")

# run 2

scope <- names(data_mag_sum %>% dplyr::select(-c(dataset.ID,RoC_max,DC,Diversity)))

compare_models(M.glm.2b,scope, order = 2)

M.glm.3a <- glm(RoC_max ~ DC * Diversity + Position, family = "gaussian" , data = data_mag_sum)
M.glm.3b <- glm(RoC_max ~ DC * Diversity * Position, family = "gaussian" , data = data_mag_sum)

anova(M.glm.2b,M.glm.3a, test= "Chisq")
anova(M.glm.3a,M.glm.3b, test= "Chisq")

# run 3

M.glm.4a <- glm(RoC_max ~ DC * Diversity * Position + SMOOTH, family = "gaussian" , data = data_mag_sum)
M.glm.4b <- glm(RoC_max ~ DC * Diversity * Position * SMOOTH, family = "gaussian" , data = data_mag_sum)


anova(M.glm.3b,M.glm.4a, test= "Chisq")
anova(M.glm.4a,M.glm.4b, test= "Chisq")

# glm.fin

M.glm.fin <- M.glm.4b

# ----------------------------------------------
#
#               FIGURES 
#
# ----------------------------------------------

Color.legen_01 <- brewer.pal(n = levels(data_mag_summmary$SMOOTH) %>% length(), name = 'Set2')
names(Color.legen_01) <- levels(data_mag_summmary$SMOOTH)


Color.legen_02 <- brewer.pal(n = levels(data_mag_summmary$DC) %>% length(), name = 'Set3')
names(Color.legen_02) <- levels(data_mag_summmary$DC)

Color.legen_03 <- brewer.pal(n = levels(data_mag_summmary$Position) %>% length(), name = 'Set1')
names(Color.legen_03) <- levels(data_mag_summmary$Position)

Color.legen_04 <- brewer.pal(n = levels(data_mag_summmary$Diversity) %>% length(), name = 'Paired')
names(Color.legen_04) <- levels(data_mag_summmary$Diversity)


##############
#   FIG 1
#############

data_example <- list(dataset.id = tibble_Europe2$dataset.id[[2]],
                     filtered.counts = tibble_Europe2$filtered.counts[[2]],
                     list_ages = tibble_Europe2$list_ages[[2]])

data_example_MW <- fc_calculate_RoC_comparison(data_example,
                                               Working.Unit = "MW",
                                               BIN.size = 500,
                                               N.shifts = 5,
                                               rand = 1000,
                                               interest.treshold =  age_lim)
data_example_MW <- within(data_example_MW, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
data_example_MW <- within(data_example_MW, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))

FIG1_visual_example_MW <- data_example_MW %>%
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
  xlab("Age (cal yr BC)")+ylab("Rate of Change score")+
  facet_grid(SMOOTH~DC, scales = "free_x")

FIG1_visual_example_MW

ggsave("~/RESULTS/Methods/FIN/FIG1_visual_example_MW.pdf",
       plot = FIG1_visual_example_MW,
       width = 20, height = 12, units = "cm")


data_example_BIN <- fc_calculate_RoC_comparison(data_example,
                                               Working.Unit = "BINs",
                                               BIN.size = 500,
                                               rand = 1000,
                                               interest.treshold =  age_lim)
data_example_BIN <- within(data_example_BIN, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
data_example_BIN <- within(data_example_BIN, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))

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
  xlab("Age (cal yr BC)")+ylab("Rate of Change score")+
  facet_grid(SMOOTH~DC, scales = "free_x")

FIG1_visual_example_BIN

ggsave("~/RESULTS/Methods/FIN/FIG1_visual_example_BIN.pdf",
       plot = FIG1_visual_example_BIN,
       width = 20, height = 12, units = "cm")



##############
#   FIG 2
#############
M.glm.fin

formula(M.glm.fin)

data_mag_summmary <- data_mag_sum %>%
  group_by(DC, Diversity, Position, SMOOTH) %>%
  summarise(N= n(),
            RoC_max.m = mean(RoC_max),
            RoC_max_SD =  sd(RoC_max, na.rm = T),
            RoC_max_SE =  sd(RoC_max)/sqrt(N),
            RoC_max_05 = quantile(RoC_max,0.025),
            RoC_max_95 = quantile(RoC_max,0.975)
            )


data_mag_summmary <- within(data_mag_summmary, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
data_mag_summmary <- within(data_mag_summmary, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))

data_mag_summmary <- within(data_mag_summmary, Position <- factor(Position, levels = c("recent","late")))
levels(data_mag_summmary$Position) <- c("high levels density","low level density")

data_mag_summmary <- within(data_mag_summmary, Diversity <- factor(Diversity, levels = c("low","high")))
levels(data_mag_summmary$Diversity) <- c("low diversity","high diversity")

FIG2_mag_MW  <- data_mag_summmary %>%
  ggplot(aes(y=RoC_max.m,x=Position, fill=SMOOTH))+
  geom_bar(stat="identity", position = "dodge", color="gray30")+
  geom_errorbar(aes(ymax=RoC_max.m+RoC_max_SE, ymin=RoC_max.m-RoC_max_SE), 
                position = position_dodge(width=0.9), width=0.2, size=0.5, color="gray50")+
  facet_grid(DC~Diversity, scales = "free_y")+
  theme_classic()+
  scale_fill_manual(values = Color.legen_01)+
  theme(legend.position = "bottom")+
  labs(x= "",
       y= "Mean maximum Rate-of-Change score",
       fill = "Smoothing")

FIG2_mag_MW

ggsave("~/RESULTS/Methods/FIN/FIG2_mag_MW.pdf",
       plot = FIG2_mag_MW,
       height = 15, width = 20, units="cm")

##############
#   FIG 3
#############

formula(glm.fin)

data_success_summary <- data_success_sum %>%
  filter(SEGMENT == "focus") %>%
  ungroup() %>%
  group_by(Position, SMOOTH, Diversity) %>%
  summarise(VALUE.M = mean(VALUE, na.rm = T),
            VALUE.SD = sd(VALUE, na.rm = T),
            VALUE.SE = sd(VALUE, na.rm = T)/sqrt(n()),
            VALUE.05 = quantile(VALUE,0.025, na.rm = T),
            VALUE.95 = quantile(VALUE,0.975, na.rm = T)
  ) %>%
  ungroup()



FIG3_sum_MW_gam_A  <-data_success_summary %>%
  ggplot(aes(y=VALUE.M,x=SMOOTH,fill=SMOOTH, group=SMOOTH))+
  geom_bar(stat="identity", position="dodge", color="gray30")+
  geom_errorbar(aes(ymin=VALUE.M-VALUE.SE,ymax=VALUE.M+VALUE.SE, group=SMOOTH),
                position=position_dodge(width=0.9), width=0.2, size=0.5, color="gray50")+
  facet_grid(Diversity~Position)+
  labs(y="Percentage of Peak detection",
       x="Smoothing",
       fill = "Smoothing",
       title = "focal area (correct detection)")+
  theme_classic()+
  coord_cartesian(ylim=c(0,0.75))+
  scale_fill_manual(values = Color.legen_01)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

FIG3_sum_MW_gam_A

formula(glm.fin.e)

data_success_summary_E <- data_success_sum %>%
  filter(SEGMENT == "empty") %>%
  ungroup() %>%
  group_by(Position, SMOOTH, Diversity, DC) %>%
  summarise(VALUE.M = mean(VALUE, na.rm = T),
            VALUE.SD = sd(VALUE, na.rm = T),
            VALUE.SE = sd(VALUE, na.rm = T)/sqrt(n()),
            VALUE.05 = quantile(VALUE,0.025, na.rm = T),
            VALUE.95 = quantile(VALUE,0.975, na.rm = T)
  ) %>%
  ungroup()


FIG3_sum_MW_gam_B  <-data_success_summary_E %>%
  ggplot(aes(y=VALUE.M,x=SMOOTH,fill=DC, group=DC))+
  geom_bar(stat="identity", position="dodge", color="gray30")+
  geom_errorbar(aes(ymin=VALUE.M-VALUE.SE,ymax=VALUE.M+VALUE.SE, group=DC),
                position=position_dodge(width=0.9), width=0.2, size=0.5, color="gray50")+
  facet_grid(Diversity~Position)+
  labs(y="Percentage of Peak detection",
       x="Smoothing",
       fill = "DC",
       title = "outside of focal area (false positive)")+
  theme_classic()+
  coord_cartesian(ylim=c(0,0.3))+
  scale_fill_manual(values = Color.legen_02)+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right")

FIG3_sum_MW_gam_B

FIG3_sum_MW_gam <- ggarrange(
  FIG3_sum_MW_gam_A,
  FIG3_sum_MW_gam_B,
  nrow = 1, widths = c(0.8,1)
)

FIG3_sum_MW_gam

FIG3_sum_MW_gam_fin <- annotate_figure(FIG3_sum_MW_gam, left = "Percentage of Peak detection", bottom = " Smoothing") 

FIG3_sum_MW_gam_fin

ggsave("~/RESULTS/Methods/FIN/FIG3_sum_MW_gam.pdf",
       plot = FIG3_sum_MW_gam_fin,
       height = 12, width = 25, units="cm")



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
  filter(data_site_A$list_ages$ages$age < 8000) %>%
  dim()

fc_get_pollen_data(data_site_A, sm.type = "shep",N.taxa = 80) %>%
  filter(age < 8000) %>%
  group_by(name) %>%
  summarise(SUM = sum(value)) %>%
  arrange(-SUM)
  

#Site B


which(tibble_Europe2$dataset.id %in%  4012 )

data_site_B <- list(dataset.id = tibble_Europe2$dataset.id[[45]],
                    filtered.counts = tibble_Europe2$filtered.counts[[45]],
                    list_ages = tibble_Europe2$list_ages[[45]])

data_site_B_pollen <-fc_get_pollen_data(data_site_B, sm.type = "shep",N.taxa = 10)

data_site_B$filtered.counts %>%
  as_tibble() %>%
  filter(data_site_B$list_ages$ages$age < 8000) %>%
  dim()

fc_get_pollen_data(data_site_B, sm.type = "shep",N.taxa = 50) %>%
  filter(age < 8000) %>%
  group_by(name) %>%
  summarise(SUM = sum(value)) %>%
  arrange(-SUM)


# SITE c

which(tibble_Europe2$dataset.id %in%  40951 )

data_site_C <- list(dataset.id = tibble_Europe2$dataset.id[[224]],
                    filtered.counts = tibble_Europe2$filtered.counts[[224]],
                    list_ages = tibble_Europe2$list_ages[[224]])


data_site_C_pollen <-fc_get_pollen_data(data_site_C, sm.type = "shep",N.taxa = 10)

data_site_C$filtered.counts %>%
  filter(data_site_C$list_ages$ages$age < 8000) %>%
  as_tibble() %>%
  dim()

fc_get_pollen_data(data_site_C, sm.type = "shep",N.taxa = 103) %>%
  filter(age < 8000) %>%
  group_by(name) %>%
  summarise(SUM = sum(value)) %>%
  arrange(-SUM)


# Site D

which(tibble_Europe2$dataset.id %in%  4314 )

data_site_D <- list(dataset.id = tibble_Europe2$dataset.id[[50]],
                    filtered.counts = tibble_Europe2$filtered.counts[[50]],
                    list_ages = tibble_Europe2$list_ages[[50]])

data_site_D$filtered.counts %>%
  filter(data_site_D$list_ages$ages$age < 8000) %>%
  rowSums()

data_site_D_pollen <-fc_get_pollen_data(data_site_D, sm.type = "shep",N.taxa = 10)

data_site_D$filtered.counts %>%
  as_tibble() %>%
  filter(data_site_D$list_ages$ages$age < 8000) %>%
  dim()

fc_get_pollen_data(data_site_D, sm.type = "shep",N.taxa = 33) %>%
  filter(age < 8000) %>%
  group_by(name) %>%
  summarise(SUM = sum(value)) %>%
  arrange(-SUM)


# Pollen taxa table

common_taxa<- c(data_site_A_pollen$name,
                data_site_B_pollen$name,
                data_site_C_pollen$name,
                data_site_D_pollen$name
                ) %>%
  unique()

library (RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
Palette.1<- getPalette(length(common_taxa))
names(Palette.1)<- sort(common_taxa)  


# Rate of Change

data_site_A_RoC <- fc_R_ratepol(data.source.pollen = data_site_A$filtered.counts,
                              data.source.age = data_site_A$list_ages,
                              sm.type = "shep", 
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              Working.Unit = "MW",
                              BIN.size = 500,
                              N.shifts = 5,
                              rand = 1000,
                              standardise = T, 
                              S.value = 150, 
                              DC = "euc.sd",
                              interest.treshold = age_lim,
                              Debug = F)


data_Site_B_RoC <- fc_R_ratepol(data.source.pollen = data_site_B$filtered.counts,
                              data.source.age = data_site_B$list_ages,
                              sm.type = "shep", 
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              Working.Unit = "MW",
                              BIN.size = 500,
                              N.shifts = 5,
                              rand = 1000,
                              standardise = T, 
                              S.value = 150, 
                              DC = "euc.sd",
                              interest.treshold = age_lim,
                              Debug = F)


data_Site_C_RoC <- fc_R_ratepol(data.source.pollen = data_site_C$filtered.counts,
                              data.source.age = data_site_C$list_ages,
                              sm.type = "shep", 
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              Working.Unit = "MW",
                              BIN.size = 500,
                              N.shifts = 5,
                              rand = 1000,
                              standardise = T, 
                              S.value = 150, 
                              DC = "euc.sd",
                              interest.treshold = age_lim,
                              Debug = F)


data_Site_D_RoC <- fc_R_ratepol(data.source.pollen = data_site_D$filtered.counts,
                                data.source.age = data_site_D$list_ages,
                              sm.type = "shep", 
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              Working.Unit = "MW",
                              BIN.size = 500,
                              N.shifts = 5,
                              rand = 1000,
                              standardise = T, 
                              S.value = 150, 
                              DC = "euc.sd",
                              interest.treshold = age_lim,
                              Debug = F)

# FIGURES 

FIG4_Site_A_1 <- data_site_A$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = c(0,3e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  geom_rug(sides = "b")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,3e-4))+
  labs(x="Age (cal yr BC)",
      y="Density of samples"
       )+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
    )+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))


FIG4_Site_B_1 <- data_site_B$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = c(0,3e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  geom_rug(sides = "b")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,3e-4))+
  labs(x="Age (cal yr BC)",
      y="Density of samples"
  )+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )+
coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))

FIG4_Site_C_1 <- data_site_C$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = c(0,3e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  geom_rug(sides = "b")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  scale_y_continuous(breaks = c(0,3e-4))+
  labs(x="Age (cal yr BC)",
       y="Density of samples"
  )+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )+
coord_flip(xlim = c(age_lim,0), ylim = c(0,3e-4))


FIG4_Site_D_1 <- data_site_D$list_ages$ages %>%
  filter(age < 9000) %>%
  ggplot(aes(x=age))+
  geom_hline(yintercept = c(0,3e-4), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_density(color="gray30", fill="gray50")+
  geom_rug(sides = "b")+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  labs(x="Age (cal yr BC)",
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


FIG4_Site_A_2 <-  data_site_A_pollen %>%
  bind_rows(.,data.frame(sample.id=NA,name=common_taxa,value=0,depth=NA, age = -100, newage=NA)) %>%
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
    x="Age (cal yr BC)",
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

FIG4_Site_B_2 <-  data_site_B_pollen %>%
  bind_rows(.,data.frame(sample.id=NA,name=common_taxa,value=0,depth=NA, age = -100, newage=NA)) %>%
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
    x="Age (cal yr BC)",
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

FIG4_Site_C_2 <-  data_site_C_pollen %>%
  bind_rows(.,data.frame(sample.id=NA,name=common_taxa,value=0,depth=NA, age = -100, newage=NA)) %>%
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
    x="Age (cal yr BC)",
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

FIG4_Site_D_2 <-  data_site_D_pollen %>%
  bind_rows(.,data.frame(sample.id=NA,name=common_taxa,value=0,depth=NA, age = -100, newage=NA)) %>%
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
    x="Age (cal yr BC)",
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

FIG4_Site_A_3 <- data_site_A_RoC %>%
ggplot(aes(y=ROC, 
           x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,20))+
  geom_hline(yintercept = seq(from=0,to=20, by=5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BC)",
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


FIG4_Site_B_3 <- data_Site_B_RoC %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,20))+
  geom_hline(yintercept = seq(from=0,to=20, by=5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BC)",
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

FIG4_Site_C_3 <- data_Site_C_RoC %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,20))+
  geom_hline(yintercept = seq(from=0,to=20, by=5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BC)",
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


FIG4_Site_D_3 <- data_Site_D_RoC %>%
  ggplot(aes(y=ROC, 
             x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  coord_flip(xlim = c(age_lim,0), ylim = c(0,20))+
  geom_hline(yintercept = seq(from=0,to=20, by=5), color="gray80", size=0.1)+
  geom_vline(xintercept = seq(from=0,to=age_lim, by=2000), color="gray80", size=0.1)+
  geom_ribbon(aes(ymin=ROC.dw, ymax=ROC.up), alpha=1/2, color="gray80", fill="gray80")+
  geom_line(alpha=1, size=0.5)+
  geom_point(data = . %>% filter(PEAK ==T ),color="green", size=2, shape=16, alpha=2/3)+
  geom_hline(yintercept = 0, color="purple", size=0.1)+
  labs(
    x="Age (cal yr BC)",
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


data_site_A_RoC %>%
  filter(PEAK == T)

data_Site_B_RoC %>%
  filter(PEAK == T)

data_Site_C_RoC %>%
  filter(PEAK == T)

data_Site_D_RoC %>%
  filter(PEAK == T)



library(cowplot)
FIG4_Site_A <- plot_grid(FIG4_Site_A_1, FIG4_Site_A_2, FIG4_Site_A_3,
                         align = "h",axis = "bt", ncol = 3, rel_widths = c(1,3,1))
FIG4_Site_B <- plot_grid(FIG4_Site_B_1, FIG4_Site_B_2, FIG4_Site_B_3,
                         align = "h",axis = "bt", ncol = 3, rel_widths = c(1,3,1))
FIG4_Site_C <- plot_grid(FIG4_Site_C_1, FIG4_Site_C_2, FIG4_Site_C_3,
                         align = "h",axis = "bt", ncol = 3, rel_widths = c(1,3,1))
FIG4_Site_D <- plot_grid(FIG4_Site_D_1, FIG4_Site_D_2, FIG4_Site_D_3,
                         align = "h",axis = "bt", ncol = 3, rel_widths = c(1,3,1))


FIG4_Site_comparison <- ggarrange(
  FIG4_Site_A,
  FIG4_Site_B,
  FIG4_Site_C,
  FIG4_Site_D,
  #labels = c(17334,"","",40951,"","",4012,"",""),
  labels = c("A","B","C","D"),
  heights = c(1.8,1,1,1),
  ncol = 1, nrow = 4, legend = "none")
FIG4_Site_comparison

FIG4_Site_comparison_legend <- ggarrange(
  FIG4_Site_comparison,
  my_legend, 
  ncol=1, heights = c(12,2))
FIG4_Site_comparison_legend

ggsave("~/RESULTS/Methods/FIN/FIG4_Site_comparison.pdf",
       plot = FIG4_Site_comparison, 
       height = 20, width = 22, units="cm")

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
                                     #axis.ticks.x = element_blank(),
                                     #axis.text.x = element_blank(),
                                     legend.position = "none"
                                     )+
                                #ylab("")+
                                ylab("value of env. variable")+
                                xlab("Age (cal yr BC)"),
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
                               ylab("% of pollen grains")+
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

Supplementary_F1


ggsave("~/RESULTS/Methods/FIN/Supplementary_F1.pdf",
       plot = Supplementary_F1,
       height = 10, width = 15, units="cm")


##########
# FIG S2 #
########## 

data_supp_fig_MW <- rbind(
  perform_sim_ld_recent_MW,
  perform_sim_ld_late_MW,
  perform_sim_hd_recent_MW,
  perform_sim_hd_late_MW
)

data_supp_fig_MW$Dataset_type <- c(rep("LD-R",120),rep("LD-L",120),rep("HD-R",120),rep("HD-L",120))

data_supp_fig_MW <- within(data_supp_fig_MW, DC <- factor(DC, levels = c("euc","euc.sd","chord","chisq")))
data_supp_fig_MW <- within(data_supp_fig_MW, SMOOTH <- factor(SMOOTH, levels = c("none","m.avg","grim","age.w","shep")))
data_supp_fig_MW <- within(data_supp_fig_MW, Dataset_type <- factor(Dataset_type, levels = c("LD-R","LD-L","HD-R","HD-L")))
data_supp_fig_MW <- within(data_supp_fig_MW, SEGMENT <- factor(SEGMENT, levels = c("focus","empty")))
data_supp_fig_MW <- within(data_supp_fig_MW, PEAK <- factor(PEAK, levels = c("PEAK.T","PEAK.G","PEAK.S")))
levels(data_supp_fig_MW$PEAK) <- c("Threshold","GAM","SNI")


Supplementary_F2 <- data_supp_fig_MW %>%
  ggplot(aes(y=VALUE.M,x=Dataset_type,fill=SEGMENT, group=SEGMENT))+
  geom_bar(stat="identity", position="dodge", color="gray30")+
  geom_errorbar(aes(ymin=VALUE.M-VALUE.SD,ymax=VALUE.M+VALUE.SD, group=SEGMENT),
                position=position_dodge(width=0.9), width=0.2, size=0.5, color="gray50")+
  facet_grid(SMOOTH~DC+PEAK)+
  scale_fill_manual("Position in sequence", labels=c("focal area (correct detection)","outside of focal area (false positive)"),
                    values = c("darkseagreen","coral"))+
  ylab("Percentage of Peak detection")+xlab("Type of simulated dataset")+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = -45,hjust = 0.3, vjust = 0.2 ))


Supplementary_F2

ggsave("~/RESULTS/Methods/FIN/Supplementary_F2.pdf",
       plot = Supplementary_F2,
       height = 15, width = 22, units="cm")

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
    scale_fill_manual(values = Color.legen_03)+
    coord_cartesian(ylim=c(0,1))+
    labs(x="desnisty of levels",
         y="Percentage of Peak detection"
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
    scale_fill_manual(values = Color.legen_01)+
    coord_cartesian(ylim=c(0,1))+
    labs(x="Smoothing",
         y="Percentage of Peak detection")+
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
    scale_fill_manual(values=Color.legen_04)+
    coord_cartesian(ylim=c(0,1))+
    labs(x= "Diversity of pollen taxa",
         y= "Percentage of Peak detection")+
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
    scale_fill_manual(values = Color.legen_02)+
    coord_cartesian(ylim=c(0,1))+
    labs(x= "Dissimilarity coeficient",
         y= "Percentage of Peak detection")+
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

Success_supp_A_a <- annotate_figure(Success_supp_A, top = "focal area (correct detection)")

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
    scale_fill_manual(values = Color.legen_03)+
    coord_cartesian(ylim=c(0,1))+
    labs(x="desnisty of levels",
         y="Percentage of Peak detection"
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
    scale_fill_manual(values = Color.legen_01)+
    coord_cartesian(ylim=c(0,1))+
    labs(x="Smoothing",
         y="Percentage of Peak detection")+
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
    scale_fill_manual(values=Color.legen_04)+
    coord_cartesian(ylim=c(0,1))+
    labs(x= "Diversity of pollen taxa",
         y= "Percentage of Peak detection")+
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
    scale_fill_manual(values = Color.legen_02)+
    coord_cartesian(ylim=c(0,1))+
    labs(x= "Dissimilarity coeficient",
         y= "Percentage of Peak detection")+
    theme(legend.position = "none", 
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  , nrow = 1, ncol = 4
)


Success_supp_B

Success_supp_B_a <- annotate_figure(Success_supp_B, top="outside of focal area (false positive)")

Success_supp_B_a

Success_supp <- ggarrange(
  Success_supp_A_a,
  Success_supp_B_a,
  nrow = 2
)

Success_supp

Supplementary_F3 <- annotate_figure(Success_supp, left = "Percentage of Peak detection")

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
  scale_fill_manual(values = Color.legen_02)+
  labs(x= "Smoothing",
       y= "Percentage of Peak detection",
       title = "outside of focal area (false positive)")+
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

####################################
#               SAVE               #
####################################

# save.image("~/DATA/temp/ENV_METHOD_20200525.RData")

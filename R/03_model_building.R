#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                     Model building
#
#                     Ondrej Mottl 
#                         2020
#
#----------------------------------------------------------#

# load config 
source("R/00_config.R")

#----------------------------------------------------------#
# 1. Load data  -----
#----------------------------------------------------------#

perform_sim <-  read_rds("data/output/datasets/success_rate/sim_success.rds") 

#----------------------------------------------------------#
# 2. Prepare data -----
#----------------------------------------------------------#

# adjust all variables to levels and save as new dataframe
data_sum <- 
  perform_sim$raw_data %>%  
  mutate(
    success = as.numeric(success),
    
    position = factor(
      position,
      levels = c("breaks_recent", "breaks_late")) %>%
      fct_recode(
        "high density level" = "breaks_recent",
        "low density level" = "breaks_late"),
    diversity = factor(
      diversity,
      levels = c("low_diversity", "high_diversity")) %>% 
      fct_recode(
        "low richness" = "low_diversity",
        "high richness" = "high_diversity"),
    smooth = factor(
      smooth,
      levels = c("none", "m.avg", "grim", "age.w", "shep")) %>% 
      fct_recode(
        "None" = "none",
        "M_avg" = "m.avg",
        "Grimm" = "grim",
        "Age_w" = "age.w",
        "Shep" = "shep"),
    DC = factor(
      DC,
      levels = c("chord","chisq")) %>%
      fct_recode(
        "Chord" = "chord",
        "Chisq" = "chisq"),
    WU = factor(
      WU,
      levels = c("levels", "bins", "MW")),
    Peak = factor(
      Peak,
      levels = c(
        "Peak_threshold",
        "Peak_trend_linear",
        "Peak_trend_non_linear",
        "Peak_trend_GAM_deriv",
        "Peak_SNI")) %>% 
      fct_recode(
        "threshold" = "Peak_threshold",
        "trend linear" = "Peak_trend_linear",
        "trend non linear" = "Peak_trend_non_linear",
        "GAM first deriv" = "Peak_trend_GAM_deriv",
        "SNI" = "Peak_SNI"),
    dataset_ID = as.factor(dataset_ID),
    calculation_ID = as.factor(calculation_ID),
    dataset_type =  paste0(position," - ",diversity) %>% 
      as.factor(),
    RoC_setting = paste0(smooth," - ",DC) %>% 
      as.factor()) 


# explore data
data_sum %>% 
  summary()

DataExplorer::plot_str(data_sum)
DataExplorer::plot_intro(data_sum)


#----------------------------------------------------------#
# 3. Model fitting  -----
#----------------------------------------------------------#

#---------------------------------------------#
# 3.1 Peak points & WU   -----
#---------------------------------------------#

#--------------------------------#
# 3.1.1 Sucess -----
#--------------------------------#

# subset data to only include succesfull detection and adjust the succes for beta family
data_correct <- 
  data_sum %>%
  filter(segment == "focus") %>%
  ungroup() %>%
  dplyr::select(-c(segment))%>% 
  mutate(success = ifelse(success == 0, success + very_small_value, success)) %>% 
  mutate(success = ifelse(success == 1, success - very_small_value, success)) 

summary(data_correct)

nrow(data_correct)

write_rds(data_correct, "data/output/datasets/subsets/data_correct.rds")

# fit actual model

start_t <- Sys.time()

mod_correct <-  
  glmmTMB(success ~ WU * Peak * dataset_type + # fixed effets 
            (1|dataset_ID) +(1|RoC_setting), # random effects
          data = data_correct,
          family = beta_family(link = "logit")
  )
end_t <- Sys.time()

end_t-start_t

summary(mod_correct)
# r2(mod_correct)
# check_distribution(mod_correct)
# check_singularity(mod_correct)
# check_heteroscedasticity(mod_correct)
# check_model(mod_correct)
# model_performance(mod_correct)
# check_autocorrelation(mod_correct)
qplot(residuals(mod_correct))

write_rds(mod_correct,"data/output/models/mod_correct.rds")

# set up cluster
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, c("data_correct", "n_cores"), envir = environment())
parallel::clusterEvalQ(cl,library("glmmTMB"))

c(getAllTerms(mod_correct))

# estime th best model by AIC
mod_correct_dd <- 
  pdredge(
    mod_correct,
    subset = "cond(dataset_type)",
    cluster = cl,
    trace = T)

mod_correct_dd %>% 
  as_tibble() %>% 
  filter(delta < 2) %>% 
  View()

mod_correct_dd %>%
  as_tibble() %>%
  write_csv("data/output/result_tables/mod_correct_compare.csv")

# -> full model
mod_correct_select <- mod_correct
write_rds(mod_correct_select,"data/output/models/mod_correct_select.rds")

#--------------------------------#
# 3.1.2 false posities  -----
#--------------------------------#

# subset data to only include false positive detection and adjust the succes for beta family
data_false <- 
  data_sum %>%
  filter(segment == "empty") %>%
  ungroup() %>%
  dplyr::select(-c(segment))%>% 
  mutate(success = ifelse(success == 0, success + very_small_value, success)) %>% 
  mutate(success = ifelse(success == 1, success - very_small_value, success)) 


summary(data_false)

nrow(data_false)

write_rds(data_false, "data/output/datasets/subsets/data_false.rds")

# fit actual model

start_t <- Sys.time()

mod_false <-  
  glmmTMB(success ~ WU * Peak * dataset_type + # fixed effets 
            (1|dataset_ID) +(1|RoC_setting), # random effects
          data = data_false,
          family = beta_family(link = "logit")
  )
end_t <- Sys.time()

end_t-start_t

summary(mod_false)
# r2(mod_false)
# check_distribution(mod_false)
# check_singularity(mod_false)
# check_heteroscedasticity(mod_false)
# check_model(mod_false)
# model_performance(mod_false)
# check_autocorrelation(mod_false)
qplot(residuals(mod_false))

write_rds(mod_false,"data/output/models/mod_false.rds")


# set up cluster
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, c("data_false", "n_cores"), envir = environment())
parallel::clusterEvalQ(cl,library("glmmTMB"))

c(getAllTerms(mod_false))

# estime th best model by AIC
mod_false_dd <- 
  pdredge(
    mod_false,
    subset = "cond(dataset_type)",
    cluster = cl,
    trace = T)


mod_false_dd %>% 
  as_tibble() %>% 
  filter(delta < 2) %>% 
  View()


mod_false_dd %>%
  as_tibble() %>%
  write_csv("data/output/result_tables/mod_false_compare.csv")

# -> full model
mod_false_select <- mod_false
write_rds(mod_false_select,"data/output/models/mod_false_select.rds")

#---------------------------------------------#
# 3.2 RoC and dataset properties   -----
#---------------------------------------------#

#--------------------------------#
# 3.2.1 Sucess -----
#--------------------------------#

# subset data to only include succesfull detection and adjust the succes for beta family
data_detail_correct <- 
  data_sum %>%
  filter(segment == "focus") %>%
  filter(Peak == "trend non linear") %>% 
  filter(WU == "MW") %>% 
  ungroup() %>%
  dplyr::select(-c(segment, Peak, WU, dataset_type, RoC_setting)) %>% 
  mutate(success = ifelse(success == 0, success + very_small_value, success)) %>% 
  mutate(success = ifelse(success == 1, success - very_small_value, success)) 


summary(data_detail_correct)

nrow(data_detail_correct)

write_rds(data_detail_correct, "data/output/datasets/subsets/data_detail_correct.rds")

# fit the model
start_t <- Sys.time()

mod_detail_correct <-  
  glmmTMB(success ~ diversity * position *  smooth *  DC + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_correct,
          family = beta_family(link = "logit")
  )
end_t <- Sys.time()

end_t-start_t

summary(mod_detail_correct)
# r2(mod_detail_correct)
# check_distribution(mod_detail_correct)
# check_singularity(mod_detail_correct)
# check_heteroscedasticity(mod_detail_correct)
# check_model(mod_detail_correct)
# model_performance(mod_detail_correct)
# check_autocorrelation(mod_detail_correct)
qplot(residuals(mod_detail_correct))


write_rds(mod_detail_correct,"data/output/models/mod_detail_correct.rds")

# set up cluster
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, c("data_detail_correct", "n_cores"), envir = environment())
parallel::clusterEvalQ(cl,library("glmmTMB"))

# estime th best model by AIC
mod_detail_correct_dd <- 
  pdredge(
    mod_detail_correct,
    cluster = cl,
    trace = T)

mod_detail_correct_dd %>%
  as_tibble() %>%
  write_csv("data/output/result_tables/mod_detail_correct_compare.csv")


# more models have similar AIC
mod_detail_correct_dd %>% 
  as_tibble() %>% 
  filter(delta < 2) %>% 
  View()

mod_detail_correct_m1 <-  
  glmmTMB(success ~ DC + position +  smooth + 
            DC:position + DC:smooth + position:smooth + 
            DC:position:smooth +  # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_correct,
          family = beta_family(link = "logit")
  )


mod_detail_correct_m2 <-  
  glmmTMB(success ~ DC + diversity + position +  smooth + 
            DC:position + DC:smooth + position:smooth + 
            DC:position:smooth + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_correct,
          family = beta_family(link = "logit")
  )

mod_detail_correct_final_comp <-
  compare_performance(
    mod_detail_correct_m1,
    mod_detail_correct_m2,
    rank = T)

mod_detail_correct_final_comp

mod_detail_correct_final_comp %>% 
  as_tibble() %>% 
  write_csv(
    .,"data/output/result_tables/mod_detail_correct_final_comp.csv")

# -> model 1 is the best
mod_detail_correct_select <-  mod_detail_correct_m1
write_rds(mod_detail_correct_select,"data/output/models/mod_detail_correct_select.rds")

#--------------------------------#
# 3.2.2 false positives -----
#--------------------------------#


# subset data to only include succesfull detection and adjust the succes for beta family
data_detail_false <- 
  data_sum %>%
  filter(segment == "empty") %>%
  filter(Peak == "trend non linear") %>% 
  filter(WU == "MW") %>% 
  ungroup() %>%
  dplyr::select( -c(segment, Peak, WU, dataset_type, RoC_setting)) %>% 
  mutate(
    success = ifelse(
      success == 0,
      success + very_small_value,
      success
      )
    ) %>% 
  mutate(
    success = ifelse(
      success == 1,
      success - very_small_value,
      success)) 


summary(data_detail_false)

nrow(data_detail_false)

write_rds(data_detail_false, "data/output/datasets/subsets/data_detail_false.rds")

# fit the model
start_t <- Sys.time()

mod_detail_false <-  
  glmmTMB(success ~ diversity * position *  smooth *  DC + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit")
  )
end_t <- Sys.time()

end_t-start_t

summary(mod_detail_false)
# r2(mod_detail_false)
# check_distribution(mod_detail_false)
# check_singularity(mod_detail_false)
# check_heteroscedasticity(mod_detail_false)
# check_model(mod_detail_false)
# model_performance(mod_detail_false)
# check_autocorrelation(mod_detail_false)
qplot(residuals(mod_detail_false))


write_rds(mod_detail_false,"data/output/models/mod_detail_false.rds")


# set up cluster
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
parallel::clusterExport(
  cl, 
  c("data_detail_false", "n_cores"),
  envir = environment())
parallel::clusterEvalQ(cl,library("glmmTMB"))


# estime th best model by AIC
mod_detail_false_dd <- 
  pdredge(
    mod_detail_false,
    cluster = cl,
    trace = T)

# more models have similar AIC
mod_detail_false_dd %>% 
  as_tibble() %>% 
  filter(delta < 2) %>% 
  View()

mod_detail_false_dd %>%
  as_tibble() %>%
  write_csv("data/output/result_tables/mod_detail_false_compare.csv")

# select the best model

mod_detail_false_m1 <- 
  glmmTMB(success ~ DC + position +  smooth + 
            DC:smooth + position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))

mod_detail_false_m2 <- 
  glmmTMB(success ~ DC + position +  smooth + 
            DC:position + DC:smooth + position:smooth + 
            DC:position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))
  

mod_detail_false_m3 <- 
  glmmTMB(success ~ DC + diversity + position +  smooth + 
            DC:diversity + DC:smooth + position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))

mod_detail_false_m4 <- 
  glmmTMB(success ~ DC + diversity + position +  smooth + 
            DC:diversity + DC:position + DC:smooth + position:smooth + 
            DC:position:smooth +
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))

mod_detail_false_m5 <- 
  glmmTMB(success ~ DC + position +  smooth + 
            DC:position + DC:smooth + position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))

mod_detail_false_m6 <- 
  glmmTMB(success ~ DC + diversity + position +  smooth + 
            DC:diversity + DC:position + DC:smooth + position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))

mod_detail_false_m7 <- 
  glmmTMB(success ~ DC + diversity + position +  smooth + 
            DC:smooth + position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))

mod_detail_false_m8 <- 
  glmmTMB(success ~ DC + diversity + position +  smooth + 
            DC:diversity + DC:smooth + diversity:position + position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))

mod_detail_false_m9 <- 
  glmmTMB(success ~  position +  smooth + 
            position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))

mod_detail_false_m10 <- 
  glmmTMB(success ~ DC + position +  smooth + 
            position:smooth + 
            (1|dataset_ID), # random effects
          data = data_detail_false,
          family = beta_family(link = "logit"))


mod_detail_false_final_comp <-
  compare_performance(
    mod_detail_false_m1,
    mod_detail_false_m2,
    mod_detail_false_m3,
    mod_detail_false_m4,
    mod_detail_false_m5,
    mod_detail_false_m6,
    mod_detail_false_m7,
    mod_detail_false_m8,
    mod_detail_false_m9,
    mod_detail_false_m10,
    rank = T)

mod_detail_false_final_comp

mod_detail_false_final_comp %>% 
  as_tibble() %>% 
  write_csv(
    .,"data/output/result_tables/mod_detail_false_final_comp.csv")




# -> model 2 is the best
mod_detail_false_select <-  mod_detail_correct_m2

write_rds(mod_detail_false_select,"data/output/models/mod_detail_false_select.rds")

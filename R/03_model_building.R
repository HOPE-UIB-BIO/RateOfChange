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

perform_sim <-  read_rds("data/output/datasets/simulated_success.rds") 

#----------------------------------------------------------#
# 2. Prepare data -----
#----------------------------------------------------------#

# adjust all variables to levels and save as new dataframe
data_success_sum <- 
  perform_sim$raw_data %>%  
  mutate(
    success = as.numeric(success),
    
    position = factor(
      position,
      levels = c("breaks_recent", "breaks_late")) %>%
      fct_recode(
        "high density level" = "breaks_recent",
        "low_density level" = "breaks_late"),
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
        "M.avg" = "m.avg",
        "Grimm" = "grim",
        "Age.w" = "age.w",
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
    dataset_type =  paste0(position,"_",diversity) %>% 
      as.factor(),
    RoC_setting = paste0(smooth,"_",DC) %>% 
      as.factor()) 

data_success_sum %>% 
  summary()

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
data_success_focus <- 
  data_success_sum %>%
  filter(segment == "focus") %>%
  ungroup() %>%
  dplyr::select(-c(segment))%>% 
  mutate(success = ifelse(success == 0, success + very_small_value, success)) %>% 
  mutate(success = ifelse(success == 1, success - very_small_value, success)) 


summary(data_success_focus)

nrow(data_success_focus)

write_rds(data_success_focus, "data/output/datasets/data_success_focus.rds")

# fit actual model

start_t <- Sys.time()

mod_success_focus_full <-  
  glmmTMB(success ~ WU * Peak * dataset_type + # fixed effets 
            (1|dataset_ID) +(1|RoC_setting), # random effects
          data = data_success_focus,
          family = beta_family(link = "logit")
  )
end_t <- Sys.time()

end_t-start_t

summary(mod_success_focus_full)
r2(mod_success_focus_full)
check_distribution(mod_success_focus_full)
check_singularity(mod_success_focus_full)
check_heteroscedasticity(mod_success_focus_full)
check_model(mod_success_focus_full)
model_performance(mod_success_focus_full)
check_autocorrelation(mod_success_focus_full)
qplot(residuals(mod_success_focus_full))

write_rds(mod_success_focus_full,"data/output/models/temp_model_success_focus.rds")

# set up cluster
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl); 
parallel::clusterExport(cl,c("data_success_focus","n_cores"),envir=environment());
parallel::clusterEvalQ(cl,library("glmmTMB"))

getAllTerms(mod_success_focus_full)

# estime th best model by AIC
mod_success_focus_dd <- 
  pdredge(
    mod_success_focus,
    subset = "cond(dataset_type)",
    cluster = cl,
    trace = T)

mod_success_focus_dd %>% 
  as_tibble() %>% 
  filter(delta < 2) %>% 
  View()

mod_success_focus_dd %>%
  as_tibble() %>%
  write_csv("data/output/result_tables/mod_success_focus.csv")

# -> full model
mod_success_focus_select <- mod_success_focus_full
write_rds(mod_success_focus_select,"data/output/models/temp_model_success_focus_select.rds")

#--------------------------------#
# 3.1.2 false posities  -----
#--------------------------------#


# subset data to only include false positive detection and adjust the succes for beta family
data_success_empty <- 
  data_success_sum %>%
  filter(segment == "empty") %>%
  ungroup() %>%
  dplyr::select(-c(segment))%>% 
  mutate(success = ifelse(success == 0, success + very_small_value, success)) %>% 
  mutate(success = ifelse(success == 1, success - very_small_value, success)) 


summary(data_success_empty)

nrow(data_success_empty)

write_rds(data_success_empty, "data/output/datasets/data_success_empty.rds")

# fit actual model

start_t <- Sys.time()

mod_success_empty_full <-  
  glmmTMB(success ~ WU * Peak * dataset_type + # fixed effets 
            (1|dataset_ID) +(1|RoC_setting), # random effects
          data = data_success_empty,
          family = beta_family(link = "logit")
  )
end_t <- Sys.time()

end_t-start_t

summary(mod_success_empty_full)
r2(mod_success_empty_full)
check_distribution(mod_success_empty_full)
check_singularity(mod_success_empty_full)
check_heteroscedasticity(mod_success_empty_full)
check_model(mod_success_empty_full)
model_performance(mod_success_empty_full)
check_autocorrelation(mod_success_empty_full)
qplot(residuals(mod_success_empty_full))

write_rds(mod_success_empty_full,"data/output/models/temp_model_success_empty.rds")


# set up cluster
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl); 
parallel::clusterExport(cl,c("data_success_empty","n_cores"),envir=environment());
parallel::clusterEvalQ(cl,library("glmmTMB"))

getAllTerms(mod_success_empty_full)

# estime th best model by AIC
mod_success_empty_dd <- 
  pdredge(
    mod_success_empty_full,
    subset = "cond(dataset_type)",
    cluster = cl,
    trace = T)


mod_success_empty_dd %>% 
  as_tibble() %>% 
  filter(delta < 2) %>% 
  View()


mod_success_empty_dd %>%
  as_tibble() %>%
  write_csv("data/output/result_tables/mod_success_empty.csv")

# -> full model
mod_success_empty_select <- mod_success_empty_full
write_rds(mod_success_empty_select,"data/output/models/temp_model_success_empty_select.rds")

#---------------------------------------------#
# 3.2 RoC and dataset properties   -----
#---------------------------------------------#

#--------------------------------#
# 3.2.1 Sucess -----
#--------------------------------#

# subset data to only include succesfull detection and adjust the succes for beta family
data_detail_focus <- 
  data_success_sum %>%
  filter(segment == "focus") %>%
  filter(Peak == "trend non linear") %>% 
  filter(WU == "MW") %>% 
  ungroup() %>%
  dplyr::select(-c(segment, Peak, WU, dataset_type, RoC_setting)) %>% 
  mutate(success = ifelse(success == 0, success + very_small_value, success)) %>% 
  mutate(success = ifelse(success == 1, success - very_small_value, success)) 


summary(data_detail_focus)

nrow(data_detail_focus)

write_rds(data_detail_focus, "data/output/datasets/data_detail_focus.rds")

# fit the model
start_t <- Sys.time()

mod_detail_focus <-  
  glmmTMB(success ~ diversity * position *  smooth *  DC + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_focus,
          family = beta_family(link = "logit")
  )
end_t <- Sys.time()

end_t-start_t

summary(mod_detail_focus)
r2(mod_detail_focus)
check_distribution(mod_detail_focus)
check_singularity(mod_detail_focus)
check_heteroscedasticity(mod_detail_focus)
check_model(mod_detail_focus)
model_performance(mod_detail_focus)
check_autocorrelation(mod_detail_focus)
qplot(residuals(mod_detail_focus))


write_rds(mod_detail_focus,"data/output/models/temp_model_detail_focus.rds")

# set up cluster
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl); 
parallel::clusterExport(cl,c("data_detail_focus","n_cores"),envir=environment());
parallel::clusterEvalQ(cl,library("glmmTMB"))

getAllTerms(data_detail_focus)

# estime th best model by AIC
mod_detail_focus_dd <- 
  pdredge(
    mod_detail_focus,
    cluster = cl,
    trace = T)

mod_detail_focus_dd %>%
  as_tibble() %>%
  write_csv("data/output/result_tables/mod_detail_focus.csv")


# more models have similar AIC
mod_detail_focus_dd %>% 
  as_tibble() %>% 
  filter(delta < 2) %>% 
  View()

mod_detail_focus_m1 <-  
  glmmTMB(success ~ position +  smooth + position:smooth + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_focus,
          family = beta_family(link = "logit")
  )


mod_detail_focus_m2 <-  
  glmmTMB(success ~ DC + position +  smooth  + position:smooth + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_focus,
          family = beta_family(link = "logit")
  )

mod_detail_focus_m3 <-  
  glmmTMB(success ~ diversity + position + smooth  + position:smooth + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_focus,
          family = beta_family(link = "logit")
  )

compare_performance(
  mod_detail_focus_m1,
  mod_detail_focus_m2,
  mod_detail_focus_m3,
  rank = T
)

# -> model 1 is the best
mod_detail_focus_select <-  mod_detail_focus_m1
write_rds(mod_detail_focus_select,"data/output/models/temp_model_detail_focus_select.rds")

#--------------------------------#
# 3.2.2 false positives -----
#--------------------------------#


# subset data to only include succesfull detection and adjust the succes for beta family
data_detail_empty <- 
  data_success_sum %>%
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


summary(data_detail_empty)

nrow(data_detail_empty)

write_rds(data_detail_empty, "data/output/datasets/data_detail_empty.rds")

# fit the model
start_t <- Sys.time()

mod_detail_empty <-  
  glmmTMB(success ~ diversity * position *  smooth *  DC + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_empty,
          family = beta_family(link = "logit")
  )
end_t <- Sys.time()

end_t-start_t

summary(mod_detail_empty)
r2(mod_detail_empty)
check_distribution(mod_detail_empty)
check_singularity(mod_detail_empty)
check_heteroscedasticity(mod_detail_empty)
check_model(mod_detail_empty)
model_performance(mod_detail_empty)
check_autocorrelation(mod_detail_empty)
qplot(residuals(mod_detail_empty))


write_rds(mod_detail_empty,"data/output/models/temp_model_detail_empty.rds")


# set up cluster
cl <- parallel::makeCluster(n_cores)

doParallel::registerDoParallel(cl)

parallel::clusterExport(
  cl, 
  c("data_detail_empty", "n_cores"),
  envir = environment())

parallel::clusterEvalQ(cl,library("glmmTMB"))

getAllTerms(data_detail_empty)

# estime th best model by AIC
mod_detail_empty_dd <- 
  pdredge(
    mod_detail_empty,
    cluster = cl,
    trace = T)

# more models have similar AIC
mod_detail_empty_dd %>% 
  as_tibble() %>% 
  filter(delta < 2) %>% 
  View()

mod_detail_empty_dd %>%
  as_tibble() %>%
  write_csv("data/output/result_tables/mod_detail_empty.csv")


mod_detail_empty_m1 <-  
  glmmTMB(success ~ position +  smooth + position:smooth + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_empty,
          family = beta_family(link = "logit")
  )


mod_detail_empty_m2 <-  
  glmmTMB(success ~ DC + position +  smooth  + position:smooth + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_empty,
          family = beta_family(link = "logit")
  )

mod_detail_empty_m3 <-  
  glmmTMB(success ~ diversity + position + smooth  + position:smooth + # fixed effets
            (1|dataset_ID), # random effects
          data = data_detail_empty,
          family = beta_family(link = "logit")
  )


compare_performance(
  mod_detail_empty_m1,
  mod_detail_empty_m2,
  mod_detail_empty_m3,
  rank = T
)

# -> model 1 is the best
mod_detail_empty_select <-  mod_detail_empty_m1
write_rds(mod_detail_empty_select,"data/output/models/temp_model_detail_empty_select.rds")

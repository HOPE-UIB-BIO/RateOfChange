#################################################################
### ------------ Estimate RoC values for all sites ---------- ###
#################################################################

# ----------------------------------------------
#                     SETUP
# ----------------------------------------------
library(tidyverse)
devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package")
library(RRatepol)

# ----------------------------------------------
#                     LOAD DATA
# ----------------------------------------------

HOPE_data_smooth <- readRDS("C:/Users/ondre/OneDrive - University of Bergen/HOPE_data/R_Data/_HOPE_DATA/HOPE_data_smooth20200910.RDS")

glimpse(HOPE_data_smooth)

Region_filter_tibble <- tibble(REGION = c("Europe","North America","Asia","Latin America","Africa","Oceania"),
                               max_limit =c(rep(8.5e3,2),rep(12e3,4)))

HOPE_data_work <- HOPE_data_smooth %>% 
  left_join(Region_filter_tibble, by="REGION")


# ----------------------------------------------
#               COMPUTATION 
# ----------------------------------------------

s.time <- Sys.time()

HOPE_data_work_ROC <-  HOPE_data_work %>%
  mutate(., ROC = pmap(list(filtered.counts,list_ages,max_limit),
                       .f = function(.x,.y,.z)
                         {
                         try(res <- RRatepol::fc_estimate_RoC(
                           data_source_pollen = .x,
                           data_source_age = .y,
                           smooth_method  = "age.w", 
                           Working_Units = "MW",
                           bin_size  = 500,
                           Number_of_shifts  = 5,
                           rand = 1000,
                           standardise = T,
                           DC = "chisq",
                           interest_threshold  = .z,
                           Peak = "GAM",
                           Debug = F))
                         
                         if(!exists("res")){res <- "NA"}
                         return(res)
                         } ))

f.time <- Sys.time()
tot.time <- f.time - s.time
tot.time


# ----------------------------------------------
#                 SAVE RESULT 
# ----------------------------------------------

HOPE_data_work_ROC <- HOPE_data_work_ROC %>%
  dplyr::select(dataset.id, ROC) %>%
  filter(purrr::map(ROC, is_tibble)==T)


saveRDS(HOPE_data_work_ROC, file = "~/HOPE/DATA/output/HOPE_Roc20200917.RDS") 

# ----------------------------------------------
#               CLEAN UP 
# ----------------------------------------------
rm(list = ls())

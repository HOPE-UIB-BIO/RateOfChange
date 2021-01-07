#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                   Data preparation
#
#                     Ondrej Mottl 
#                         2020
#
#----------------------------------------------------------#

# load config 
source("R/00_config.R")

#----------------------------------------------------------#
# 1. Simulate datasets -----
#----------------------------------------------------------#

# low diversity recent

sim_ld_recent <- 
  .simulate.pollen.data.in.multiple.datasets(
    time = time_seq, 
    nforc = N_env, 
    nprox = low_diversity, 
    manual_edit = TRUE,
    breaks = breaks_recent,
    jitter = TRUE,
    rarity = TRUE,
    N_datasets = N_rep)


sim_ld_late <- 
  .simulate.pollen.data.in.multiple.datasets(
    time=time_seq, 
    nforc=N_env, 
    nprox=high_diversity, 
    manual_edit = T,
    breaks=breaks_late,
    jitter = T,
    rarity=T,
    N_datasets=N_rep)


sim_hd_recent <- 
  .simulate.pollen.data.in.multiple.datasets(
    time=time_seq, 
    nforc=N_env, 
    nprox=high_diversity, 
    manual_edit = T,
    breaks=breaks_recent,
    jitter = T,
    rarity=T,
    N_datasets=N_rep)


sim_hd_late <- 
  .simulate.pollen.data.in.multiple.datasets(
    time=time_seq, 
    nforc=N_env, 
    nprox=high_diversity, 
    manual_edit = T,
    breaks=breaks_late,
    jitter = T,
    rarity=T,
    N_datasets=N_rep)


#----------------------------------------------------------#
# 2. Merge datasets -----
#----------------------------------------------------------#

simulated_dataset <- 
  dpyr::bind_rows(
    
    tibble::tibble(
      sim_ld_recent,
      diversity = "low_diversity",
      position = "breaks_recent"),
    
    tibble::tibble(
      sim_ld_late,
      diversity = "low_diversity",
      position = "breaks_late"),
    
    tibble::tibble(
      sim_ld_recent,
      diversity = "high_diversity",
      position = "breaks_recent"),
    
    tibble::tibble(
      sim_ld_late,
      diversity = "high_diversity",
      position = "breaks_late")) %>% 
  mutate(
    dataset_ID = paste(ID,diversity, position) %>% 
      as.factor() %>% 
      as.numeric()) %>% 
  arrange(dataset_ID) %>% 
  dplyr::select(dataset_ID, diversity, position, community_data, list_ages) %>% 
  mutate(dataset_ID = as.character(dataset_ID)) 

simulated_dataset

#----------------------------------------------------------#
# 3. Save datasets -----
#----------------------------------------------------------#
write_rds(
  simulated_dataset,
  "data/output/datasets/simulated/simulated_dataset.rds")

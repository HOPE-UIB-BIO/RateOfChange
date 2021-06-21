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

# low diversity early

sim_ld_early <- 
  .simulate.pollen.data.in.multiple.datasets(
    time = time_seq, 
    nforc = N_env, 
    nprox = low_diversity, 
    manual_edit = TRUE,
    breaks = breaks_early,
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


sim_hd_early <- 
  .simulate.pollen.data.in.multiple.datasets(
    time=time_seq, 
    nforc=N_env, 
    nprox=high_diversity, 
    manual_edit = T,
    breaks=breaks_early,
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
      sim_ld_early,
      diversity = "low_diversity",
      position = "breaks_early"),
    
    tibble::tibble(
      sim_ld_late,
      diversity = "low_diversity",
      position = "breaks_late"),
    
    tibble::tibble(
      sim_ld_early,
      diversity = "high_diversity",
      position = "breaks_early"),
    
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
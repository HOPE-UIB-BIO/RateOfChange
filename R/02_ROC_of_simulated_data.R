#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                 RoC in simulated datasets
#                 & peak detection success
#
#                     Ondrej Mottl 
#                         2020
#
#----------------------------------------------------------#

# load config 
source("R/00_config.R")

#----------------------------------------------------------#
# 1. Load data -----
#----------------------------------------------------------#

list_files_output <-  list.files("data/output/datasets/")

if(any(list_files_output %in% "simulated_dataset.rds")){
  simulated_dataset <-  read_rds("data/output/datasets/simulated_dataset.rds") 
} else {
  source("R/01_data_creation.R")
}

#----------------------------------------------------------#
# 2. Calculate RoC -----
#----------------------------------------------------------#

sim_ROC_levels <- 
  fc_estimate_RoC_by_all_methods(
    simulated_dataset,
    Working_Unit = "levels", 
    interest_threshold = 8000)

write_rds(
  sim_ROC_levels,
  "data/output/datasets/sim_ROC_levels_compress.rds",
  compress = "gz")

sim_ROC_bins <- 
  fc_estimate_RoC_by_all_methods(
    simulated_dataset,
    Working_Unit = "bins",
    bin_size = 500, 
    Number_of_shifts = 1,
    interest_threshold = 8000)

write_rds(
  sim_ROC_bins,
  "data/output/datasets/sim_ROC_bins_compress.rds",
  compress = "gz")


sim_ROC_MW <- 
  fc_estimate_RoC_by_all_methods(
    simulated_dataset,
    Working_Unit = "MW",
    bin_size = 500, 
    Number_of_shifts = 5,
    interest_threshold = 8000)

write_rds(
  sim_ROC_MW,
  "data/output/datasets/sim_ROC_MW_compress.rds",
  compress = "gz")


#----------------------------------------------------------#
# 3. Merge files -----
#----------------------------------------------------------#
sim_ROC_all <-
  bind_rows(
    tibble(sim_ROC_levels, WU = "levels"),
    tibble(sim_ROC_bins, WU = "bins"),
    tibble(sim_ROC_MW, WU = "MW")
  ) %>% 
  mutate(
    calculation_ID = paste0(WU, calculation_number) %>% 
      as.factor() %>% 
      as.numeric())  %>% 
  dplyr::select(dataset_ID, calculation_ID, everything()) %>% 
  arrange(dataset_ID, calculation_ID)

sim_ROC_all$calculation_ID %>% 
  unique() %>% 
  length()

sim_ROC_all$dataset_ID %>% 
  unique() %>% 
  length()

write_rds(
  sim_ROC_all,
  "data/output/datasets/simulated_roc.rds",
  compress = "gz")

#----------------------------------------------------------#
# 4. Detect sucesss of peak points  -----
#----------------------------------------------------------#

perform_sim <-  fc_test_success_in_simulated_data(sim_ROC_all)

write_rds(perform_sim, "data/output/datasets/simulated_success.rds")


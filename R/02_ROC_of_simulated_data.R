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

list_files_output <-  list.files("data/output/datasets/simulated/")

if(any(list_files_output %in% "simulated_dataset.rds")){
  simulated_dataset <-  
    read_rds("data/output/datasets/simulated/simulated_dataset.rds") 
} else {
  source("R/01_data_creation.R")
}

#----------------------------------------------------------#
# 2. Calculate RoC -----
#----------------------------------------------------------#

ROC_levels <- 
  .estimate.RoC.by.all.methods(
    simulated_dataset,
    Working_Unit = "levels", 
    interest_threshold = 8000)

write_rds(
  ROC_levels,
  "data/output/datasets/RoC/ROC_levels_compress.rds",
  compress = "gz")

ROC_bins <- 
  .estimate.RoC.by.all.methods(
    simulated_dataset,
    Working_Unit = "bins",
    bin_size = 500, 
    Number_of_shifts = 1,
    interest_threshold = 8000)

write_rds(
  ROC_bins,
  "data/output/datasets/RoC/ROC_bins_compress.rds",
  compress = "gz")


ROC_MW <- 
  .estimate.RoC.by.all.methods(
    simulated_dataset,
    Working_Unit = "MW",
    bin_size = 500, 
    Number_of_shifts = 5,
    interest_threshold = 8000)

write_rds(
  ROC_MW,
  "data/output/datasets/RoC/ROC_MW_compress.rds",
  compress = "gz")


#----------------------------------------------------------#
# 3. Merge files -----
#----------------------------------------------------------#
ROC_all <-
  bind_rows(
    tibble(ROC_levels, WU = "levels"),
    tibble(ROC_bins, WU = "bins"),
    tibble(ROC_MW, WU = "MW")
  ) %>% 
  mutate(
    calculation_ID = paste0(WU, calculation_number) %>% 
      as.factor() %>% 
      as.numeric())  %>% 
  dplyr::select(dataset_ID, calculation_ID, everything()) %>% 
  arrange(dataset_ID, calculation_ID)

ROC_all$calculation_ID %>% 
  unique() %>% 
  length()

ROC_all$dataset_ID %>% 
  unique() %>% 
  length()

write_rds(
  ROC_all,
  "data/output/datasets/RoC/ROC_all.rds",
  compress = "gz")

#----------------------------------------------------------#
# 4. Detect sucesss of peak points  -----
#----------------------------------------------------------#

perform_sim <-  .test.success.in.simulated.data(ROC_all)

write_rds(perform_sim, "data/output/datasets/success_rate/sim_success.rds")


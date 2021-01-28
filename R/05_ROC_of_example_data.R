#----------------------------------------------------------#
#
#
#             Rate-of-change in palaeoecology 
#
#                     Calculation RoC 
#                     in example data
#
#                     Ondrej Mottl 
#                         2020
#
#----------------------------------------------------------#

# load config 
source("R/00_config.R")

#----------------------------------------------------------#
# 1. Load data and prepare data -----
#----------------------------------------------------------#

example_data <- RRatepol::example_data

data_site_A <- 
  list(dataset_ID = example_data$dataset.id[[4]],
       community_data = example_data$filtered.counts[[4]],
       age_data = example_data$list_ages[[4]]$ages,
       uncertainity_data = example_data$list_ages[[4]]$age_position)

write_rds(
  data_site_A,
  "data/output/datasets/pollen_sites/data_site_A.rds"
)

data_site_B <- 
  list(dataset_ID = example_data$dataset.id[[1]],
       community_data = example_data$filtered.counts[[1]],
       age_data = example_data$list_ages[[1]]$ages,
       uncertainity_data = example_data$list_ages[[1]]$age_position)

write_rds(
  data_site_B,
  "data/output/datasets/pollen_sites/data_site_B.rds"
)

data_site_C <- 
  list(dataset_ID = example_data$dataset.id[[2]],
       community_data = example_data$filtered.counts[[2]],
       age_data = example_data$list_ages[[2]]$ages,
       uncertainity_data = example_data$list_ages[[2]]$age_position)

write_rds(
  data_site_C,
  "data/output/datasets/pollen_sites/data_site_C.rds"
)

data_site_D <- 
  list(dataset_ID = example_data$dataset.id[[3]],
       community_data = example_data$filtered.counts[[3]],
       age_data = example_data$list_ages[[3]]$ages,
       uncertainity_data = example_data$list_ages[[3]]$age_position)

write_rds(
  data_site_D,
  "data/output/datasets/pollen_sites/data_site_D.rds"
)


# check size
data_site_A$community_data[ ,-1] %>%
  as_tibble() %>%
  dim()

data_site_B$community_data[ ,-1] %>%
  as_tibble() %>%
  dim()

data_site_C$community_data[ ,-1] %>%
  as_tibble() %>%
  dim()

data_site_D$community_data[ ,-1] %>%
  as_tibble() %>%
  dim()

#----------------------------------------------------------#
# 2. Calculate RoC  -----
#----------------------------------------------------------#

# 2.1. Site A ----- 
data_site_A_RoC_levels <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_A$community_data,
    data_source_age = data_site_A$age_data,
    age_uncertainty = data_site_A$uncertainity_data,
    smooth_method  = "shep",
    Working_Units  = "levels",
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_A_RoC_levels %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_A_RoC_levels.rds")

data_site_A_RoC_bins <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_A$community_data,
    data_source_age = data_site_A$age_data,
    age_uncertainty = data_site_A$uncertainity_data,
    smooth_method  = "shep", 
    Working_Units  = "bins",
    bin_size  = 500,
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_A_RoC_bins %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_A_RoC_bins.rds")


data_site_A_RoC_MW <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_A$community_data,
    data_source_age = data_site_A$age_data,
    age_uncertainty = data_site_A$uncertainity_data,
    smooth_method  = "shep", 
    Working_Units  = "MW",
    bin_size  = 500,
    Number_of_shifts  = 5,
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_A_RoC_MW %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_A_RoC_MW.rds")

# 2.2. Site B ----- 
data_site_B_RoC_levels <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_B$community_data,
    data_source_age = data_site_B$age_data,
    age_uncertainty = data_site_B$uncertainity_data,
    smooth_method  = "shep",
    Working_Units  = "levels",
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_B_RoC_levels %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_B_RoC_levels.rds")

data_site_B_RoC_bins <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_B$community_data,
    data_source_age = data_site_B$age_data,
    age_uncertainty = data_site_B$uncertainity_data,
    smooth_method  = "shep", 
    Working_Units  = "bins",
    bin_size  = 500,
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_B_RoC_bins %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_B_RoC_bins.rds")


data_site_B_RoC_MW <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_B$community_data,
    data_source_age = data_site_B$age_data,
    age_uncertainty = data_site_B$uncertainity_data,
    smooth_method  = "shep", 
    Working_Units  = "MW",
    bin_size  = 500,
    Number_of_shifts  = 5,
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_B_RoC_MW %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_B_RoC_MW.rds")

# 2.3. Site C ----- 
data_site_C_RoC_levels <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_C$community_data,
    data_source_age = data_site_C$age_data,
    age_uncertainty = data_site_C$uncertainity_data,
    smooth_method  = "shep",
    Working_Units  = "levels",
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_C_RoC_levels %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_C_RoC_levels.rds")

data_site_C_RoC_bins <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_C$community_data,
    data_source_age = data_site_C$age_data,
    age_uncertainty = data_site_C$uncertainity_data,
    smooth_method  = "shep", 
    Working_Units  = "bins",
    bin_size  = 500,
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_C_RoC_bins %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_C_RoC_bins.rds")


data_site_C_RoC_MW <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_C$community_data,
    data_source_age = data_site_C$age_data,
    age_uncertainty = data_site_C$uncertainity_data,
    smooth_method  = "shep", 
    Working_Units  = "MW",
    bin_size  = 500,
    Number_of_shifts  = 5,
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_C_RoC_MW %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_C_RoC_MW.rds")

# 2.4. Site D ----- 
data_site_D_RoC_levels <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_D$community_data,
    data_source_age = data_site_D$age_data,
    age_uncertainty = data_site_D$uncertainity_data,
    smooth_method  = "shep",
    Working_Units  = "levels",
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_D_RoC_levels %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_D_RoC_levels.rds")

data_site_D_RoC_bins <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_D$community_data,
    data_source_age = data_site_D$age_data,
    age_uncertainty = data_site_D$uncertainity_data,
    smooth_method  = "shep", 
    Working_Units  = "bins",
    bin_size  = 500,
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_D_RoC_bins %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_D_RoC_bins.rds")


data_site_D_RoC_MW <- 
  RRatepol::fc_estimate_RoC(
    data_source_community = data_site_D$community_data,
    data_source_age = data_site_D$age_data,
    age_uncertainty = data_site_D$uncertainity_data,
    smooth_method  = "shep", 
    Working_Units  = "MW",
    bin_size  = 500,
    Number_of_shifts  = 5,
    rand = roc_n_rand,
    standardise = T, 
    N_individuals  = pollen_grains, 
    DC = "chisq",
    treads = T,
    interest_threshold  = age_lim,
    Debug = F) %>% 
  RRatepol::fc_detect_peak_points(.,method = "trend_non_linear")

data_site_D_RoC_MW %>% 
  write_rds(
    ., "data/output/example_data_roc/data_site_D_RoC_MW.rds")

#----------------------------------------------------------#
# 3. Calculate RoC for a site A in all setting -----
#----------------------------------------------------------#

ROC_site_A_all_settings <- .calculate.roc.in.all.settings(data_site_A)

write_rds(
  ROC_site_A_all_settings,
  "data/output/datasets/RoC/ROC_site_A_all_settings.rds"
)


#----------------------------------------------------------#
# 4. Calculate RoC for a site A in various bin size -----
#----------------------------------------------------------#

bin_sizes <- seq(from = 100, to = 1e3, by = 100)

for(i in bin_sizes){
  
  df_w <- 
    RRatepol::fc_estimate_RoC(
      data_source_community = data_site_A$community_data,
      data_source_age = data_site_A$age_data,
      age_uncertainty = data_site_A$uncertainity_data,
      smooth_method  = "shep", 
      Working_Units  = "MW",
      bin_size  = i,
      Number_of_shifts  = 5,
      rand = roc_n_rand,
      standardise = T, 
      N_individuals  = pollen_grains, 
      DC = "chisq",
      treads = T,
      interest_threshold  = age_lim,
      Debug = F)
  
  if(exists("bin_size_result") == FALSE){
    bin_size_result <- tibble(df_w, bin = i)
  } else {
    bin_size_result <- 
      bind_rows(
        bin_size_result,
        tibble(df_w, bin = i))
  }
}


write_rds(
  bin_size_result,
  "data/output/datasets/RoC/bin_size_result.rds"
)



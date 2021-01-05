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


data_site_B <- 
  list(dataset_ID = example_data$dataset.id[[1]],
       community_data = example_data$filtered.counts[[1]],
       age_data = example_data$list_ages[[1]]$ages,
       uncertainity_data = example_data$list_ages[[1]]$age_position)


data_site_C <- 
  list(dataset_ID = example_data$dataset.id[[2]],
       community_data = example_data$filtered.counts[[2]],
       age_data = example_data$list_ages[[2]]$ages,
       uncertainity_data = example_data$list_ages[[2]]$age_position)


data_site_D <- 
  list(dataset_ID = example_data$dataset.id[[3]],
       community_data = example_data$filtered.counts[[3]],
       age_data = example_data$list_ages[[3]]$ages,
       uncertainity_data = example_data$list_ages[[3]]$age_position)


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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_A$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_A$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_A$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_B$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_B$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_B$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_C$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_C$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_C$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_D$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_D$community_data,
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
  RRatepol::fc_estimate_RoC(data_source_community = data_site_D$community_data,
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
# 3. Extract pollen data  -----
#----------------------------------------------------------#

# common taxa in each dataset
data_site_A_dom <- .extract.dominant.pollen.taxa(data_site_A)
data_site_B_dom <- .extract.dominant.pollen.taxa(data_site_B)
data_site_C_dom <- .extract.dominant.pollen.taxa(data_site_C)
data_site_D_dom <- .extract.dominant.pollen.taxa(data_site_D)

# Common pollen taxa in all datasets
common_taxa <- 
  c(data_site_A_dom,
    data_site_B_dom,
    data_site_C_dom,
    data_site_D_dom) %>%
  unique() %>%
  sub("/",".",.) %>%
  sub("-",".",.) %>%
  sub(")",".",.) %>%
  sub(".\\(","..",.) 

# prepare pollen data for the common taxa
data_site_A_pollen <- .get.pollen.data(data_site_A, common_taxa)
data_site_B_pollen <- .get.pollen.data(data_site_B, common_taxa)
data_site_C_pollen <- .get.pollen.data(data_site_C, common_taxa)
data_site_D_pollen <- .get.pollen.data(data_site_D, common_taxa)

#----------------------------------------------------------#
# 4. Create individual plots  -----
#----------------------------------------------------------#

# 4.1. Sample density -----
plot_site_A_density <- .plot.density.of.samples(data_site_A)
plot_site_B_density <- .plot.density.of.samples(data_site_B)
plot_site_C_density <- .plot.density.of.samples(data_site_C)
plot_site_D_density <- .plot.density.of.samples(data_site_D)

# 4.2. pollen distribution -----
plot_site_A_pollen <- 
  .plot.pollen.data(data_site_A_pollen, common_taxa, strip_names = T)

plot_site_B_pollen <- 
  .plot.pollen.data(data_site_B_pollen, common_taxa)

plot_site_C_pollen <- 
  .plot.pollen.data(data_site_C_pollen, common_taxa)

plot_site_D_pollen <- 
  .plot.pollen.data(data_site_D_pollen, common_taxa)



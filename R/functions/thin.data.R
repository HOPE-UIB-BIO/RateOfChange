.thin.data <- function(WU, size_of_bin = 1000){
  
  resolution_full <-
    RRatepol::fc_estimate_RoC(
      data_source_community = data_site_A$community_data, 
      data_source_age = data_site_A$age_data,  
      age_uncertainty = FALSE, 
      Working_Units = WU,
      bin_size = size_of_bin,
      smooth_method = "age.w", 
      DC = "chisq") %>% 
    RRatepol::fc_detect_peak_points(
      data_source = .,
      method = "trend_non_linear") %>% 
    dplyr::mutate(resolution = "full")
  
  # helper functions
  # reduces resolution of a pollen record.
  .thinner <- function(data, resolution = 100, offset = resolution / 2){
    n <- nrow(data)
    #get starting sample within offset of top of record
    id <- which.min(abs(runif(1, data$age [1], data$age [1] + offset) - data$age ))
    while(data$age [id[length(id)]] < (data$age [n] - offset)){
      #find nedatat sample closest to current sample + resolution
      id.new <- which.min(abs(data$age [id[length(id)]] + resolution - data$age ))
      # check nedatat sample is at least offset away from current
      if(data$age [id.new] < data$age [id[length(id)]] + offset){
        id.new <- id.new + 1
      }
      id <- c(id, id.new)
    }
    id <- id[id < n]
    data[id, ]  
  }
  
  resolution_thin <- 
    res %>% 
    set_names() %>% 
    map(~.thinner(data_site_A_merged, resolution = .x, offset = min(.x/2, 50))) %>% 
    map(~{
      
      age_sub <- .x %>% 
        select(sample.id, depth, age)
      
      pollen_sub <- .x %>% 
        select(-c(depth, age))
      
      RRatepol::fc_estimate_RoC(
        data_source_community = pollen_sub, 
        data_source_age = age_sub,  
        age_uncertainty = FALSE, 
        Working_Units = WU,
        bin_size = size_of_bin,
        smooth_method = "age.w", 
        DC = "chisq")
    }) %>% 
    map_dfr(~RRatepol::fc_detect_peak_points(
      data_source = .x,
      method = "trend_non_linear"),
      .id = "resolution")
  
  res <- 
    bind_rows(
      resolution_full,
      resolution_thin) %>% 
    dplyr::mutate(
      resolution = factor(
        resolution,
        levels = c("full", res)))
  
  return(res)
  
}
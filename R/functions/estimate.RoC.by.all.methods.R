.estimate.RoC.by.all.methods <- function(random_data,
                                           Working_Unit = "levels",
                                           bin_size = 500, 
                                           Number_of_shifts = 5, 
                                           interest_threshold = 8000){
  
  
  performance_smooth <-  c(rep("none", 2), rep("shep", 2), rep("m.avg", 2),rep("age.w", 2), rep("grim", 2))
  performance_DC <- c(rep(c("chord", "chisq"),5))
  N_datasets <- nrow(random_data)
  
  # unique number to track down the datasets
  calculation_number <- 1
  
  for(i in 1:N_datasets) {
    
    cat(paste0("dataset ",i,"/",N_datasets), "\n")
    
    for(j in 1:length(performance_smooth)){
      
      repeat {
        try( data_temp <- 
               fc_estimate_RoC( 
                 data_source_community =  random_data$community_data[[i]] %>% 
                   dplyr::mutate(
                     sample.id = as.character(dplyr::row_number())) %>% 
                   dplyr::relocate(sample.id),
                 data_source_age = random_data$list_ages[[i]] %>% 
                   dplyr::mutate(
                     sample.id = as.character(sample.id)),
                 age_uncertainty = F,
                 smooth_method  = performance_smooth[j],
                 smooth_N_points  = 5,
                 smooth_age_range  = 500, 
                 smooth_N_max  = 9,
                 Working_Units  = Working_Unit,
                 bin_size  = bin_size ,
                 Number_of_shifts = Number_of_shifts,
                 rand = 100,
                 treads = T,
                 standardise  = T,
                 N_individuals = 150,
                 tranform_to_proportions = T,
                 DC  = performance_DC[j],
                 interest_threshold  = interest_threshold, 
                 time_standardisation = 500,
                 Debug  = F) %>%
               as_tibble(),
             silent = T )
        if (exists("data_temp")==T) break
      }
  
      
      tibble_fin <- 
        tibble(
          dataset_ID = i,
          calculation_number = calculation_number,
          smooth = performance_smooth[j],
          DC = performance_DC[j],
          diversity = random_data$diversity[i],
          position = random_data$position[i],
          data_temp) %>% 
        mutate( # PEAK detection 
          Peak_threshold = fc_detect_peak_points(data_temp, method = "threshold") %>% 
            dplyr::select(Peak) %>% 
            pluck(1),
          
          Peak_trend_linear = fc_detect_peak_points(data_temp, method = "trend_linear") %>% 
            dplyr::select(Peak) %>% 
            pluck(1),
          
          Peak_trend_non_linear = fc_detect_peak_points(data_temp, method = "trend_non_linear") %>% 
            dplyr::select(Peak) %>% 
            pluck(1),
          
          Peak_trend_GAM_deriv = fc_detect_peak_points(data_temp, method = "GAM_deriv") %>% 
            dplyr::select(Peak) %>% 
            pluck(1),
          
          Peak_SNI = fc_detect_peak_points(data_temp, method = "SNI") %>% 
            dplyr::select(Peak) %>% 
            pluck(1),
        )
      
      
      rm (data_temp)
      
      calculation_number <-  calculation_number + 1
      
      if (i == 1 & j == 1){
        tibble_fin_sum <- tibble_fin
      } else {
        tibble_fin_sum <- rbind(tibble_fin_sum,tibble_fin)
      }
    }
  }
  return(tibble_fin_sum)
}

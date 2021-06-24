.remove.top.core <- function(WU, size_of_bin = 1000){
  
  res <- 
    top_core_vector %>% 
    set_names() %>% 
    map(~filter(data_site_A_merged, age >= .x)) %>% 
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
      .id = "top_core_age")
  
  
  
  return(res)
  
}

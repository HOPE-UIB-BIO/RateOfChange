.calculate.roc.in.all.settings <- function(data){
  
  smooth_method <- c("none", "shep", "m.avg", "age.w" ,"grim")
  DC_method <- c("chord", "chisq") 
  
  res_tibble <-
    expand.grid(smooth_method, DC_method) %>% 
    as_tibble() %>% 
    rename(
      "smooth_type" = Var1,
      "DC" = Var2) %>% 
    mutate_if(is.factor,as.character)
  
  res_tibble <-
    res_tibble %>% 
    mutate(
      ROC = purrr::map2(smooth_type, DC, function(x, y){
        
        roc_score <-
          RRatepol::fc_estimate_RoC(data_source_community = data$community_data,
                                    data_source_age = data$age_data,
                                    age_uncertainty = data$uncertainity_data,
                                    smooth_method  = x,
                                    smooth_N_points = 5,
                                    smooth_N_max = 9,
                                    smooth_age_range = 500,
                                    Working_Units  = "MW",
                                    rand = roc_n_rand,
                                    standardise = T, 
                                    N_individuals  = pollen_grains, 
                                    DC = y,
                                    treads = T,
                                    interest_threshold  = age_lim,
                                    Debug = F)
        roc_peak <-
          RRatepol::fc_detect_peak_points(roc_score,method = "trend_non_linear")
        
        return(roc_peak)
      })
    )
  
  return(res_tibble)
}


fc_test_success_in_simulated_data <- function(sim_data){
  
  unique_calculation_ID <- sim_data$calculation_ID %>%
    unique() %>%
    length()
  
  pb <- txtProgressBar(min = 0, max = unique_calculation_ID, style = 3)
  
  for (i in 1: unique_calculation_ID){
    
    setTxtProgressBar(pb, i)
    
    print(paste0(i,"/",unique_calculation_ID," ",round((i/unique_calculation_ID)*100,2),"%"))
    
    dataset_work <-
      sim_data %>% 
      filter(calculation_ID == i)
    
    # breaks definition
    break_type <-  dataset_work$position %>% unique()
    breaks <-  get(break_type)
    breaks_seq <-  c("empty",rep(c("focus","empty"),length(breaks)))
    
    
    # Peak definition
    peak_type_names <-  names(dataset_work)[grepl("Peak", names(dataset_work))]
    
    # prelocate space for results
    res_mat_temp_detect <- 
      array(
        data = NA,
        dim = c(length(breaks_seq),
                length(peak_type_names)),
        dimnames = list(breaks_seq, peak_type_names) ) 
    
    res_mat_temp_samples <- 
      tibble(
        segment = breaks_seq,
        N_samples = NA
      ) 
    
    #look in arround im the distace of the differences betwen ecollogical breaks
    window_size <- 500 #diff(breaks)[1]/2
    
    # extract values for each PEAK significance test type
    for(l in 1:length(peak_type_names)){
      
      for (k in 1:length(breaks_seq)){
        
        if(breaks_seq[k] == "empty"){
          
          if(k == 1){
            target <-  c(min(dataset_work$Age), breaks[k] - window_size)  
          } else {
            if(k == length(breaks_seq)){
              target <-  c(breaks[(k - 1)/2] + window_size, max(dataset_work$Age)) 
            } else {
              target <- c(breaks[(k - 1)/2] + window_size, breaks[(k + 1)/2] - window_size)  
            }
          }
        } else {
          target <- c(breaks[k/2] - window_size, breaks[k/2] + window_size)
        }
        
        # sort by values
        target <-sort(target)
        
        # test if there is a peak in selected signif value
        res_mat_temp_detect[k,l] <-  
          dataset_work %>%
          filter(Age >= target[1] & Age <= target[2]) %>%
          select(peak_type_names[l]) %>%
          pluck(1) %>%
          any(.,na.rm=T)
        
        
        suppressWarnings(
          res_mat_temp_samples$N_samples[k] <- 
            dataset_work %>%
            filter(Age >= target[1] & Age <= target[2]) %>%
            nrow())
        
      }
    }
    
    # save result from single randomization
    sum_temp <- 
      res_mat_temp_detect %>%
      as_tibble() %>%
      mutate(segment = breaks_seq) %>%
      mutate(N_samples = res_mat_temp_samples$N_samples) %>% 
      pivot_longer(-c(segment,N_samples)) %>%
      rename(Peak = name) %>%
      group_by(segment, Peak) %>%
      summarise(success = any(value), .groups = "keep") %>%
      mutate(
        calculation_ID = i,
        dataset_ID = unique(dataset_work$dataset_ID),
        smooth = unique(dataset_work$smooth),
        DC = unique(dataset_work$DC),
        diversity = unique(dataset_work$diversity),
        position = unique(dataset_work$position),
        WU = unique(dataset_work$WU)) %>%
      dplyr::select(dataset_ID, calculation_ID, WU, diversity, position, smooth, DC, Peak, segment, success) %>% 
      ungroup()
    
    if (i == 1){
      sum_temp_res <- sum_temp 
    } else {
      sum_temp_res <- rbind(sum_temp_res,sum_temp)
    }
    
  }
  
  close(pb)
  
  # summary of randomisation
  fin_data <- 
    sum_temp_res %>%
    group_by(WU, diversity, position, smooth, DC, Peak, segment) %>%
    summarise(
      .groups = "keep",
      success_mean = mean(success, na.rm = T)
      ) %>%
    ungroup()
  
  return(list(raw_data = sum_temp_res, sum_data = fin_data))    
}

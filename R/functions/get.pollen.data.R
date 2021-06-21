.get.pollen.data <- function (data, common_list){
  # remove the sample ID
  # if (is.numeric(unlist(data$community_data[,1]))==F){
  #   data$community_data <- data$community_data[,-1]
  # }
  
  data_ext <- RRatepol:::fc_extract_data(data$community_data,
                                        data$age_data) %>%
    RRatepol:::fc_smooth_community_data(.,smooth_method  = "none") %>%
    RRatepol:::fc_check_data(.,proportion = T)
  
  plot_data <- 
    data_ext@Community %>%
    select(any_of(common_list)) %>%
    rownames_to_column() %>%
    pivot_longer(cols = -c(rowname)) %>%
    rename(sample.id = rowname) %>%
    inner_join(.,data_ext@Age, by= "sample.id")
  
  return(plot_data)
}

fc_test_simlutated_data_magnitude <- function(sim.data) {
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  
  N.datasets <- sim.data$dataset.ID %>%
    unique() %>%
    length()
  
  pb <- txtProgressBar(min = 0, max = 20, style = 3)
  
  for(i in 1:20) {
    
    setTxtProgressBar(pb, i)
    
    temp.tibble <- tibble(dataset.ID = 1:N.datasets,SMOOTH= performance.smooth[i],DC=performance.DC[i], RoC_max = NA,RoC_upq = NA, RoC_median = NA  )
      
    for(j in 1:N.datasets){
      res.temp<- sim.data %>%
        filter(SMOOTH == performance.smooth[i] &
               DC == performance.DC[i] &
               dataset.ID == j ) %>%
        dplyr::select(ROC) %>%
        pluck(1)
      
      temp.tibble$RoC_max[j] = max(res.temp)
      temp.tibble$RoC_upq[j] = quantile(res.temp, 0.95)
      temp.tibble$RoC_median[j] = median(res.temp)
    }
    
    if(i == 1){
      res.tibble <- temp.tibble
    } else {
      res.tibble <- rbind(res.tibble, temp.tibble)
    }
  }
  
  close(pb)
  return(res.tibble)
}

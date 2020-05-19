fc_test_simlutated_data_magnitude <- function(sim.data) {
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  
  N.datasets <- sim.data$dataset.ID %>%
    unique() %>%
    length()
  
  pb <- txtProgressBar(min = 0, max = 20, style = 3)
  
  for(i in 1:20) {
    
    setTxtProgressBar(pb, i)
    
    res.vect <- vector(mode = "double",length = N.datasets)
      
    for(j in 1:N.datasets){
      res.vect[j]<- sim.data %>%
        filter(SMOOTH == performance.smooth[i] &
               DC == performance.DC[i] &
               dataset.ID == j ) %>%
        dplyr::select(ROC) %>%
        max()
    }
    
    temp.tibble <- tibble(dataset.ID = 1:N.datasets,RoC_max = res.vect,SMOOTH= performance.smooth[i],DC=performance.DC[i]  )
    
    if(i == 1){
      res.tibble <- temp.tibble
    } else {
      res.tibble <- rbind(res.tibble, temp.tibble)
    }
  }
  
  close(pb)
  return(res.tibble)
}
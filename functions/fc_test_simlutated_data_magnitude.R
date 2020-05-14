fc_test_simlutated_data_magnitude <- function(sim.data) {
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  performance.tibble <- tibble(SMOOTH=NA, DC=NA, RoC_max=NA,RoC_max_SD=NA,RoC_max_05=NA,RoC_max_95=NA, .rows=20);
  performance.tibble$SMOOTH <- performance.smooth
  performance.tibble$DC <- performance.DC
  
  N.datasets <- sim.data$dataset.ID %>%
    unique() %>%
    length()
  
  res.vect <- vector(mode = "double",length = N.datasets)
  
  for(i in 1:20) {
    
    for(j in 1:N.datasets){
      res.vect[j]<- sim.data %>%
        filter(SMOOTH == performance.smooth[i] &
               DC == performance.DC[i] &
               dataset.ID == j ) %>%
        dplyr::select(ROC) %>%
        max()
    }
    
    performance.tibble$RoC_max[i] <- mean(res.vect)
    performance.tibble$RoC_max_SD[i] <- sd(res.vect)
    performance.tibble$RoC_max_05[i] <- quantile(res.vect,0.025)
    performance.tibble$RoC_max_95[i] <- quantile(res.vect,0.975)
  }
  return(performance.tibble)
}
.simulate.pollen.data.in.multiple.datasets <- function(time=0:10e3, 
                                                       nforc=4, 
                                                       nprox=10, 
                                                       manual_edit = T,
                                                       breaks=c(2000,3000),
                                                       jitter = T,
                                                       rarity=T,
                                                       N_datasets=100){
  
  
  for(i in 1:N_datasets) {
    
    # create random data
    random_data <- 
      .simulate.pollen.data(
        time = time,
        nforc = nforc, 
        mean = 100, 
        sdev = .15,
        nprox = nprox,
        var = 20,
        range = 20,
        manual_edit = manual_edit,
        breaks = breaks,
        jitter = jitter,
        rarity = rarity,
        transform_to_counts = T,
        N_pollen_grains = 300)
    
    
    
    r_tibble <-
      tibble(
        ID = i,
        community_data = list(random_data$community_data),
        list_ages = list(random_data$list_ages))
    
    if(i == 1){
      res <-  r_tibble
    } else {
      res <-  bind_rows(res,r_tibble)
    }
    
  }
  
  return(res)
}

fc_simulate_pollen_data_in_multiple_datasets <- function(time=0:10e3, 
                                                         nforc=4, 
                                                         mean=100, 
                                                         sdev=.15, 
                                                         nprox=10, 
                                                         var=20, 
                                                         range=15,
                                                         manual.edit = T,
                                                         breaks=c(2000,3000),
                                                         jitter = T,
                                                         rarity=T,
                                                         N.datasets=100){
  
  
  for(i in 1:N.datasets) {
    
    # create random data
    random.data <- fc_simulate_pollen_data(time=time,
                                           nforc = nforc, 
                                           mean = mean, 
                                           sdev=sdev,
                                           nprox = nprox,
                                           var=var,
                                           range = range,
                                           manual.edit = manual.edit,
                                           breaks = breaks,
                                           jitter = jitter,
                                           rarity = rarity)
    
   
    
    r_tibble = tibble(ID=i,filtered.counts = list(random.data$filtered.counts), list_ages =list( random.data$list_ages))
  
    if(i == 1){
      res <- r_tibble
    } else {
      res <- bind_rows(res,r_tibble)
    }
    
  }
  
  return(res)
}

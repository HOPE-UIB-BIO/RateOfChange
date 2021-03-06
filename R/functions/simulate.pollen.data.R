.simulate.pollen.data <- function(time=5e3:0, 
                                  nforc=4, 
                                  mean=100, 
                                  sdev=.15, 
                                  nprox=10, 
                                  var=20, 
                                  range=15,
                                  manual_edit = T,
                                  breaks=c(2000,3000),
                                  transform_to_counts = T,
                                  N_pollen_grains = 300,
                                  rarity = T,
                                  jitter = T)
{
  # create enviromental variables
  forcing <-  array(0, dim = c(length(time), nforc))
  forcing[1, ] <-  rnorm(nforc, mean, sdev)
  for(i in 2:length(time)){
    forcing[i, ] <-  rnorm(nforc, forcing[i-1, ], sdev)
  }
  
  # manual edit of env.var.
  if (manual_edit == T){
    
    if(length(breaks)%%2!=0)
      stop("number of breaks must be an even number (2-6)")
    if(length(breaks)>6)
      stop("Number of breaks must be maximum of 6 (2-6)")
    
    for(l in 1:( length(breaks) - 1 ) )
    {
      if(l%%2 == 1) # odd
      {
        forcing[time > breaks[l] & time<breaks[l + 1], ] <-  forcing[time > breaks[l] & time < breaks[l + 1], ] * (1 + sdev)  
      }
      
      if(l%%2 == 0) # even
      {
        forcing[time > breaks[l] & time < breaks[l + 1], ] <-  forcing[time > breaks[l] & time < breaks[l + 1], ] * (1 - sdev) 
      }
      
    }
    
  }
  
  # smooth env.var
  forcing <-  apply(forcing, 2, FUN = function(x) {
    low <-  lowess(x, f = .05, iter = 100)
    return(low$y)
  }) 
  
  # choose random optima and ranges for the biota
  ecology <- c()
  
  for(i in 1:nprox){
    ecology[[i]] <-  list(mean = rnorm(nforc, mean, var), sd = rgamma(nforc, range, 1))
  }
  
  # reactions of the biota to the environmental changes
  proxies <-  array(1, dim = c(length(time), nprox))
  o <-  c()
  
  for(i in 1:nprox){
    
    for(j in 1:nforc){
      proxies[,i] <- proxies[,i] * dnorm(forcing[,j], ecology[[i]]$mean[j], ecology[[i]]$sd[j])
      #o[i] <- weighted.mean(time, proxies[,i])
    }
  }
  
  # order taxa by abundance
  o <-  order(colSums(proxies), decreasing = TRUE)
  proxies <-  proxies[ ,o]
  
  # decrease the abundances of rare taxa
  if(rarity==T){
    
    for(i in 1:ncol(proxies)) {
      proxies[ ,i] <-  (proxies[ ,i] / max(1, runif(1, min = i-1, max = i)) )
    }
  }
  
  # switch to proportions
  proxies <- 
    apply(proxies, 1, FUN = function(x){(x / sum(x))}) %>% 
    t()
  
  # tranform to pollen grains 
  if(transform_to_counts == T){
    proxies <-  round(proxies * N_pollen_grains)
  }
  
  # jitter the resul the pollen data
  if(jitter == T){
    proxies <-  apply(proxies, 2, FUN = function(x) jitter(x, factor = 1.5, amount = 0))
    proxies <- round(proxies)
    proxies[proxies < 0] <-  0
  }
  
  # return 
  data_source_age <-  data.frame(sample.id = c(1:length(time)), age = time)
  data_source_pollen <-  as.data.frame(proxies)
  row.names(data_source_pollen) <-  c(1:length(time))
  
  return(list(
    community_data = data_source_pollen,
    list_ages = data_source_age,
    env_var = forcing))
}
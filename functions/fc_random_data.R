fc_random_data <- function(time=5e3:0, 
                           nforc=4, 
                           mean=100, 
                           sdev=.15, 
                           nprox=10, 
                           var=20, 
                           range=15,
                           manual.edit = T,
                           breaks=c(0,2000, 4000,6000),
                           jitter = T)
{
  # create enviromental variables
  forcing <- array(0, dim=c(length(time), nforc))
  forcing[1,] <- rnorm(nforc, mean, sdev)
  for(i in 2:length(time))
    forcing[i,] <- rnorm(nforc, forcing[i-1,], sdev)
  
  
  
  
  # manual edit of env.var.
  if (manual.edit==T)
  {
    
    if(length(breaks)%%2!=0)
      stop("number of breaks must be an even number (2-6)")
    if(length(breaks)>6)
      stop("Number of breaks must be maximum of 6 (2-6)")
    
    if(length(breaks)==2){
      forcing[time>breaks[1] & time<breaks[2],] <-  forcing[time>breaks[1] & time<breaks[2],] * 1.2
    }
    
    if(length(breaks)==4) {
      forcing[time>breaks[1] & time<breaks[2],] <-  forcing[time>breaks[1] & time<breaks[2],] * 1.2
      forcing[time>breaks[2] & time<breaks[3],] <-  forcing[time>breaks[2] & time<breaks[3],] * 0.8
      forcing[time>breaks[3] & time<breaks[4],] <-  forcing[time>breaks[3] & time<breaks[4],] * 1.2  
    }
    
    if(length(breaks)==6) {
      forcing[time>breaks[1] & time<breaks[2],] <-  forcing[time>breaks[1] & time<breaks[2],] * 1.2
      forcing[time>breaks[2] & time<breaks[3],] <-  forcing[time>breaks[2] & time<breaks[3],] * 0.8
      forcing[time>breaks[3] & time<breaks[4],] <-  forcing[time>breaks[3] & time<breaks[4],] * 1.2
      forcing[time>breaks[4] & time<breaks[5],] <-  forcing[time>breaks[4] & time<breaks[5],] * 0.8
      forcing[time>breaks[5] & time<breaks[6],] <-  forcing[time>breaks[5] & time<breaks[6],] * 1.2
    }
    
    
    
  }
  
  # smooth env.var
  forcing<- apply(forcing,2, FUN = function(x) {
    low <- lowess(x,f=.25,iter=0)
    return(low$y)
  }) 
  
  # choose random optima and ranges for the biota
  ecology <- c()
  for(i in 1:nprox)
    ecology[[i]] <- list(mean=rnorm(nforc, mean, var), sd=rgamma(nforc, range, 1))
  
  # reactions of the biota to the environmental changes
  proxies <- array(1, dim=c(length(time), nprox))
  o <- c()
  for(i in 1:nprox)
  {
    for(j in 1:nforc)
      proxies[,i] <- proxies[,i] * dnorm(forcing[,j], ecology[[i]]$mean[j], ecology[[i]]$sd[j])
    o[i] <- weighted.mean(time, proxies[,i])
  }
  o <- order(o, decreasing=TRUE)
  proxies <- proxies[,o]
  
  
  # jitter the resul the pollen data
  if(jitter==T)
  {
    proxies<- apply(proxies,2,FUN= function(x) jitter(x,factor = 1.5, amount = 0))
    proxies[proxies < 0] <- 0
  }
  
  # return 
  data.source.age<- list(ages=data.frame(sample.id =c(1:length(time)), age=time),
                         age_position= matrix(time, nrow = 1))
  data.source.pollen <- as.data.frame(proxies)
  row.names(data.source.pollen) <- c(1:length(time))
  
  return(list(filtered.counts = data.source.pollen, list_ages = data.source.age))
}
  

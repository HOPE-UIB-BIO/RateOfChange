fc_smooth <- function(data.source, 
                      sm.type="none",
                      N.points = 3,
                      grim.N.min = 3, 
                      grim.N.max = 9,
                      grim.age.max = 300)
{
  # imput variables:
  # data.source - data prepared by the function of fn_extract
  # sm.type = type of smoothing applied smooting 
  #     "none"    = data will not be smoothed 
  #     "m.avg"   = moving average
  #     "grim"    = Grimm smoothing
  #     "age.w"   = age weithed 
  #
  # N.points = Number of points for moving average, need to be an odd number
  #
  # grim.N.min = minimal number of samples to look in Grimm smoothing
  # grim.N.max = maximal number of samples to look in Grimm smoothing
  # grim.N.age = minial age range for Grimm smoothing
  #
  
  
  # split data into 2 datasets
  p.counts <-  data.source$Pollen
  age <- data.source$Age   
  
  
  # ----------------------------------------------
  #               NONE SMOOTHING 
  # ----------------------------------------------
  
  if(sm.type=="none")
  {
    print("data will not be smoothed")
    return(list(Pollen=p.counts, Age=age))
  }
  
  
  # ----------------------------------------------
  #           MOVING AVERAGE SMOOTHING 
  # ----------------------------------------------
  
  
  if(sm.type=="m.avg")
  {
    # check if N.points is and odd number
    if(N.points%%2 ==0)
      stop("N.points has to be an odd number")
    
    print(paste("data will be smoothed by moving average over",N.points,"points"))
    
    N.offset <- floor(N.points/2)
    N.first <- N.offset+1
    N.last <- nrow(p.counts)-N.offset
    
    for(j in 1:ncol(p.counts)) # for every species
    {
      col.work<- p.counts[[j]] # select the species
      col.res <- rep(0,length(col.work)) # create empty vector of same lengt for values to be saved
      
      for(i in N.first:N.last)
      {
        F.low <- i-N.offset # min position to look for averaging in each step
        F.high <- i+N.offset # max position to look for averaging in each step
        col.res[i]<-mean(col.work[F.low:F.high])
      }
      
      p.counts[N.first:N.last,j]<-col.res[N.first:N.last]
      
    }
    return(list(Pollen=p.counts, Age=age))
  }
  
  
  # ----------------------------------------------
  #               GRIMm SMOOTHING 
  # ----------------------------------------------
  
  if(sm.type == "grim")
  {
    # check if grim.N.min and grim.N.max are an odd numbers
    if(grim.N.min%%2 ==0)
      stop("grim.N.min has to be an odd number")
    if(grim.N.max%%2 ==0)
      stop("grim.N.max has to be an odd number")
    
    if(grim.N.min>grim.N.max)
      stop("grim.N.max has to be biger than grim.N.min")
    
    print(paste("data will be smoothed by Grimm method with min samples",grim.N.min,
                "max samples",grim.N.max,"and max age range of",grim.age.max))
    
    
    # halve the number of points rounded down
    N.grim.N.min <- floor(grim.N.min/2) 
    N.grim.N.max <- floor(grim.N.max/2)
    
    N.first <- N.grim.N.min+1 # first posible value to look
    N.last <- nrow(p.counts)-N.grim.N.min # last posible values to look
    
    for(j in 1:ncol(p.counts)) # for every species
    {
      col.work<- p.counts[[j]] # select the species
      col.res <- rep(0,length(col.work)) # create empty vector of same lengt for values to be saved
      
      for(i in N.first:N.last) # for each point between min and max
      {
        F.low <- i-N.grim.N.min # min position to look for averaging in each step
        F.high <- i+N.grim.N.min # max position to look for averaging in each step
        
        N.active <- N.grim.N.min # length of the seach parameter (set as half of the minimal samples in teh begining)
        
        # start for serch parametr 1 and continue until distance between
        # min sample size and max sample size (both halved)
        for (k in 1:(N.grim.N.max-N.grim.N.min))  
        {
          # create new search parameter that higher by 1
          N.active.test <- N.active+1
          
          # test if this increase does not invalidate rules.
          # 1) seach parameter cannot go outside of the sample size (up or down)
          # 2) seach parameter cannot be biger than selected maximum sample sizes  
          # 3) the age diference between samples selected by the seach paramated cannot be higher than 
          #  defined max age range
          # if all of those ARE TRUE then increase the real search parameter
          
          if( i-N.active.test > 0 & 
              i+N.active.test < nrow(p.counts) & 
              N.active.test<N.grim.N.max)
            { if (abs(age$newage[i-N.active.test]-age$newage[i-N.active.test])<grim.age.max)
              {N.active <- N.active+1}
            }
          
          # adjust the points by the new seach parameters
          F.low <- i-N.active  
          F.high <- i+N.active
          
        }
        # save mean of all points in seach parameter
        col.res[i]<-mean(col.work[F.low:F.high])
        
      }
      
      # update polen values
      p.counts[N.first:N.last,j]<-col.res[N.first:N.last]
      
    }
    return(list(Pollen=p.counts, Age=age))
  }
  
  # ----------------------------------------------
  #           AGE WEIGHTED SMOOTHING 
  # ----------------------------------------------
  
  if(sm.type=="age.w")
  {
    # check if N.points is and odd number
    if(N.points%%2 ==0)
      stop("N.points has to be an odd number")
    
    print(paste("data will be smoothed by moving average over",N.points,"points"))
    
    N.offset <- floor(N.points/2)
    N.first <- N.offset+1
    N.last <- nrow(p.counts)-N.offset
    
    for(j in 1:ncol(p.counts)) # for every species
    {
      col.work<- p.counts[[j]] # select the species
      col.res <- rep(0,length(col.work)) # create empty vector of same lengt for values to be saved
      
      for(i in N.first:N.last)
      {
        F.low <- i-N.offset # min position to look for averaging in each step
        F.high <- i+N.offset # max position to look for averaging in each step
        col.res[i]<-mean(col.work[F.low:F.high])
      }
      
      p.counts[N.first:N.last,j]<-col.res[N.first:N.last]
      
    }
    return(list(Pollen=p.counts, Age=age))
  }
  
  

  # ----------------------------------------------
  #                     FIN 
  # ----------------------------------------------
  
}
fc_smooth <- function(data.source, 
                      sm.type="none",
                      N.points = 3, 
                      grim.N.max = 9,
                      range.age.max = 300,
                      Debug = F)
{
  # imput variables:
  # data.source - data prepared by the function of fn_extract
  # sm.type = type of smoothing applied smooting 
  #     "none"    = data will not be smoothed 
  #     "m.avg"   = moving average
  #     "grim"    = Grimm smoothing
  #     "age.w"   = age weithed 
  #     "shep"    = Shepard's 5-term filter
  #
  # N.points = Number of points for (need to be an odd number). Used for moving average, Grimm and Age-Weighted
  # grim.N.max = maximal number of samples to look in Grimm smoothing
  # range.age.max = maximal age range for both Grimm and Age-weight smoothing
  #
  
  # split data into 2 datasets
  p.counts <-  data.source$Pollen
  age <- data.source$Age   
  
  # check if N.points is and odd number
  if(N.points%%2 ==0)
    stop("N.points has to be an odd number")
  
  # ----------------------------------------------
  #               NONE SMOOTHING 
  # ----------------------------------------------
  
  if(sm.type=="none")
  {
    if (Debug==T){print("data will not be smoothed")}
    
    return(list(Pollen=p.counts, Age=age))
  }
  
  
  # ----------------------------------------------
  #           MOVING AVERAGE SMOOTHING 
  # ----------------------------------------------
  
  if(sm.type=="m.avg")
  {
    if(Debug==T){print(paste("data will be smoothed by moving average over",N.points,"points"))}
    
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
      
      p.counts[,j]<-col.res
      
    }
    p.counts.small <- p.counts[N.first:N.last,]
    age.small <-age[N.first:N.last,] 
    return(list(Pollen=p.counts.small, Age=age.small))
  }
  
  
  # ----------------------------------------------
  #               GRIMm SMOOTHING 
  # ----------------------------------------------
  
  if(sm.type == "grim")
  {
    # check if grim.N.max is an odd numbers
    if(grim.N.max%%2 ==0)
      stop("grim.N.max has to be an odd number")
    
    # Check if miminal number of points in not gigger than maximum
    if(N.points>grim.N.max)
      stop("grim.N.max has to be biger than N.points")
    
    if (Debug==T)
      {
      print(paste("data will be smoothed by Grimm method with min samples",N.points,
                  "max samples",grim.N.max,"and max age range of",range.age.max))
      }
        
    # halve the number of points rounded down
    N.N.points <- floor(N.points/2) 
    N.grim.N.max <- floor(grim.N.max/2)
    
    N.first <- N.N.points+1 # first posible value to look
    N.last <- nrow(p.counts)-N.N.points # last posible values to look
    
    for(j in 1:ncol(p.counts)) # for every species
    {
      col.work<- p.counts[[j]] # select the species
      col.res <- rep(0,length(col.work)) # create empty vector of same lengt for values to be saved
      
      for(i in N.first:N.last) # for each point between min and max
      {
        F.low <- i-N.N.points # min position to look for averaging in each step
        F.high <- i+N.N.points # max position to look for averaging in each step
        
        N.active <- N.N.points # length of the seach parameter (set as half of the minimal samples in teh begining)
        
        # start for serch parametr 1 and continue until distance between
        # min sample size and max sample size (both halved)
        for (k in 1:(N.grim.N.max-N.N.points))  
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
            { if (abs(age$newage[i-N.active.test]-age$newage[i-N.active.test])<range.age.max)
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
      p.counts[,j]<-col.res[]
      
    }
    p.counts.small <- p.counts[N.first:N.last,]
    age.small <-age[N.first:N.last,] 
    return(list(Pollen=p.counts.small, Age=age.small))
  }
  
  # ----------------------------------------------
  #           AGE-WEIGHTED SMOOTHING 
  # ----------------------------------------------
  
  if(sm.type=="age.w")
  {
    if(Debug==T){print(paste("data will be smoothed by age-weighed average over",N.points,"points"))}
    
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
        
        # create small df with values around observed sample (in range of offset)
        df.work <-  data.frame(values= col.work[F.low:F.high], 
                              age = age$newage[F.low:F.high], 
                              Weight=rep(1,N.points))
        
        for (k in 1:nrow(df.work))
        {
          F.age <- df.work$age[k] # age value for selected sample
          F.age.sample <- age$newage[i] # age value of obserced sample
          F.age.dist <- abs(F.age-F.age.sample) # distance between those ages
  
          # Weith of points is calculated as range.age.max / distance bewtween oldest and youngest points.
          # If cannot be smaller than 1. Values very far away from the point 
          df.work$Weight[k] <- min(c(range.age.max/F.age.dist,1))   
        }
        
        col.res[i]<-weighted.mean(df.work$values,df.work$Weight)
      }
      
      p.counts[,j]<-col.res
      
    }
    p.counts.small <- p.counts[N.first:N.last,]
    age.small <-age[N.first:N.last,] 
    return(list(Pollen=p.counts.small, Age=age.small))
  }
  
  
  # ----------------------------------------------
  #             Shepard's 5-term filter 
  # ----------------------------------------------
  
  if(sm.type=="shep")
  {
    if(Debug==T){print(paste("data will be smoothed by Shepard's 5-term filter"))}
    
    N.points <- 5
    N.offset <- floor(N.points/2)
    N.first <- N.offset+1
    N.last <- nrow(p.counts)-N.offset
    
    for(j in 1:ncol(p.counts)) # for every species
    {
      col.work <- p.counts[[j]] # select the species
      col.res <- rep(0,length(col.work)) # create empty vector of same lengt for values to be saved
      
      for(i in N.first:N.last)
      {
        # calculate the Shepard number by equasion
        w.value  <- (17.0*col.work[i] + 12*(col.work[i+1]+col.work[i-1]) - 3.0*(col.work[i+2]+col.work[i-2])) / 35.0
      # there is posibility that the resut will be smaller than zero, if that is the case use 0 instead
      if(w.value<0){w.value<-0}
      col.res[i] <- w.value
      }
      
      p.counts[,j]<-col.res
      
    }
    p.counts.small <- p.counts[N.first:N.last,]
    age.small <-age[N.first:N.last,] 
    return(list(Pollen=p.counts.small, Age=age.small))
  }
  
}
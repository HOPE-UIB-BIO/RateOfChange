fc_smooth <- function(data.source, 
                      sm.type="none",
                      m.avg.N = 3,
                      grim.N.min = 3, 
                      grim.N.max = 9,
                      grim.age.max = 150)
{
  # imput variables:
  # data.source - data prepared by the function of fn_extract
  # sm.type = type of smoothing applied smooting 
  #     "none"    = data will not be smoothed 
  #     "m.avg"   = moving average
  #     "grim"    = Grimm smoothing
  #     "age.w"   = age weithed 
  #
  # m.avg.N = Number of points for moving average, need to be an odd number
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
    # check if m.avg.N is and odd number
    if(m.avg.N%%2 ==0)
      stop("m.avg.N has to be an odd number")
    
    print(paste("data will be smoothed by moving average over",m.avg.N,"points"))
    
    N.offset <- floor(m.avg.N/2)
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
      
      p.counts[[j]]<-col.res
      
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
    
    
    N.grim.N.min <- floor(grim.N.min/2)
    N.grim.N.max <- floor(grim.N.max/2)
    N.grim.age.max <- floor(grim.age.max/2)
  
    N.first <- N.grim.N.min+1
    N.last <- nrow(p.counts)-N.grim.N.min
    
    for(j in 1:ncol(p.counts)) # for every species
    {
      col.work<- p.counts[[j]] # select the species
      col.res <- rep(0,length(col.work)) # create empty vector of same lengt for values to be saved
      
      for(i in N.first:N.last)
      {
        F.low <- i-N.grim.N.min # min position to look for averaging in each step
        F.high <- i+N.grim.N.min # max position to look for averaging in each step
        
        N.active <- N.grim.N.min # length of the seach parameter
        
        for (k in 1:(N.grim.age.max-N.grim.N.min))
        {
            
          if( i-N.active > 0 & 
              i+N.active < nrow(p.counts) & 
              N.active<N.grim.N.max)
            { if (abs(age$newage[i-N.active]-age$newage[i-N.active])<grim.age.max)
              {N.active <- N.active+1}
            }
          
          F.low <- i-N.active
          F.high <- i+N.active
          
        }
        
        col.res[i]<-mean(col.work[F.low:F.high])
        
      }
      
      p.counts[[j]]<-col.res
      
    }
    return(list(Pollen=p.counts, Age=age))
  }
  
  # ----------------------------------------------
  #               X SMOOTHING 
  # ----------------------------------------------
  
  
}
fc_smooth <- function(data.source, 
                      sm.type="none",
                      N.points = 3, 
                      grim.N.max = 9,
                      range.age.max = 300,
                      Round.result = T,
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
  
  
  # ----------------------------------------------
  #                     SETUP 
  # ----------------------------------------------
  # split data into 2 datasets
  p.counts <-  as.data.frame(data.source$Pollen)
  age <- as.data.frame(data.source$Age)   
  focus.par <- matrix(data=NA,nrow=nrow(age),ncol=2)
  
  
  
  # ----------------------------------------------
  #               NONE SMOOTHING 
  # ----------------------------------------------
  
  if(sm.type=="none")
  {
    return(list(Pollen=p.counts, Age=age, Age.un=data.source$Age.un, Dim.val= data.source$Dim.val))
  }
  
  # ----------------------------------------------
  #                 SMOOTHING 
  # ----------------------------------------------
  
  # check if N.points is and odd number
  if(N.points%%2 ==0)
    stop("N.points has to be an odd number")
  
  # check if grim.N.max is an odd numbers
  if(grim.N.max%%2 ==0)
    stop("grim.N.max has to be an odd number")
  
  # Check if miminal number of points in not gigger than maximum
  if(N.points>grim.N.max)
    stop("grim.N.max has to be biger than N.points")
  
  
  if (Debug==T & sm.type == "none")
    {cat("data will not be smoothed",fill=T)}
  if(Debug==T & sm.type == "m.avg")
    {cat(paste("data will be smoothed by moving average over",N.points,"points"),fill=T)}
  if (Debug==T & sm.type == "grim")
    {print(cat("data will be smoothed by Grimm method with min samples",N.points,
              "max samples",grim.N.max,"and max age range of",range.age.max),fill=T)}
  if(Debug==T & sm.type == "age.w")
  {cat(paste("data will be smoothed by age-weighed average over",N.points,"points"),fill=T)}
  if(Debug==T & sm.type == "shep"){cat(paste("data will be smoothed by Shepard's 5-term filter"),fill=T)}
  
  # crete support fucntion for GRIMM smoothing
  search.parameter <- function(A, B, range.age.max)
    {
      # test if this increase does not invalidate rules.
      # 1) seach parameter cannot go outside of the sample size (up or down)
      # 2) seach parameter cannot be biger than selected maximum sample sizes  
      # 3) the age diference between samples selected by the seach paramated cannot be higher than 
      #  defined max age range
      # if all of those ARE TRUE then increase the real search parameter
      
      
      for (k in 1:(grim.N.max-N.points))  
      {
        # create new search parameter that is lower by 1
        A.test <- A-1
        if( A.test > 0 &  B-A.test < grim.N.max) # i+N.active.test < nrow(p.counts) &
        { if (abs(age$newage[A.test]-age$newage[B])<range.age.max)
        {A <- A.test}
        }
        
        # create new search parameter that higher by 1
        B.test <- B+1
        if( B.test < nrow(p.counts) &  B-A.test < grim.N.max) 
        { if (abs(age$newage[A]-age$newage[B.test])<range.age.max)
        {B <- B.test}
        }
      }
      return(c(A,B))
    }
  
 # ----------------------------------------------
 #                   CALCULATION 
 # ----------------------------------------------
  
  for(j in 1:ncol(p.counts)) # for every species
  {
    col.work <- .subset2(p.counts,j) # select the species
    col.res <- rep(0,length(col.work)) # create empty vector of same lengt for values to be saved
    
    for(i in 1:nrow(p.counts)) # for each sample
    {
      
      # ----------------------------------------------
      #           MOVING AVERAGE SMOOTHING 
      # ----------------------------------------------
      if(sm.type=="m.avg")
      {
        if( i < round(0.5*(N.points))+1 ) {  # Samples near beginning (moving window truncated)
          focus.par[i,] = c(1, ( i + round(0.5*(N.points)) ))
        } else {
          if( i > nrow(age)-round(0.5*(N.points)) ) { # Samples near end
            focus.par[i,] = c( (i - round(0.5*(N.points))), nrow(age) )  
          } else { 
            focus.par[i,] = c( (i - round(0.5*(N.points))), (i+round(0.5*(N.points))) )
          }
        }
        col.res[i]<- mean (col.work[focus.par[i,1]:focus.par[i,2]])
      }
      
      
      # ----------------------------------------------
      #               GRIMm SMOOTHING 
      # ----------------------------------------------
      if(sm.type == "grim")
      {
        if( i < round(0.5*(grim.N.max))+1 ) {  # Samples near beginning (moving window truncated)
          focus.par[i,1] <- 1 
          focus.par[i,2] <- ( i + round(0.5*(N.points)) )
          focus.par[i,] <- search.parameter(focus.par[i,1] ,focus.par[i,2], range.age.max  )
        } else {
          if( i > nrow(age)-round(0.5*(N.points)) ) { # Samples near end
            focus.par[i,1] <- (i - round(0.5*(N.points)))
            focus.par[i,2] <- nrow(age) 
            focus.par[i,] <- search.parameter(focus.par[i,1] ,focus.par[i,2], range.age.max  )
            
          } else { 
            focus.par[i,1] <- (i - round(0.5*(N.points))) 
            focus.par[i,2] <- (i + round(0.5*(N.points)))
            focus.par[i,] <- search.parameter(focus.par[i,1] ,focus.par[i,2], range.age.max  )
          }
        }
        col.res[i]<- mean(col.work[focus.par[i,1]:focus.par[i,2]])
      }
      
      # ----------------------------------------------
      #           AGE-WEIGHTED SMOOTHING 
      # ----------------------------------------------
      if(sm.type=="age.w")
      {
        
        if( i < round(0.5*(N.points))+1 ) {  # Samples near beginning (moving window truncated)
          focus.par[i,] = c(1, ( i + round(0.5*(N.points)) ))
        } else {
          if( i > nrow(age)-round(0.5*(N.points)) ) { # Samples near end
            focus.par[i,] = c( (i - round(0.5*(N.points))), nrow(age) )  
          } else { 
            focus.par[i,] = c( (i - round(0.5*(N.points))), (i+round(0.5*(N.points))) )
          }
        }
        
      # create small df with values around observed sample (in range of offset)
      df.work <-  data.frame(values= col.work[focus.par[i,1]:focus.par[i,2]],
                             age = age$newage[focus.par[i,1]:focus.par[i,2]], 
                             Weight=1)
      
      # Weith of points is calculated as range.age.max / distance bewtween oldest and youngest points.
      # If cannot be smaller than 1. Values very far away from the point 
      F.age.dist <- abs(df.work$age-age$newage[i])
      const <-  range.age.max/F.age.dist
      const[const>1] <- 1
      df.work$Weight <- const
      
      col.res[i]<-weighted.mean(df.work$values,df.work$Weight)
      
      }
      
      # ----------------------------------------------
      #             Shepard's 5-term filter 
      # ----------------------------------------------
      if(sm.type=="shep")
      {
        if(i < round(0.5*(N.points))+1) {
          col.res[i] <- col.work[i]
        } else {
          if (i > nrow(age)-round(0.5*(N.points)) ) {
            col.res[i] <- col.work[i]
          } else {
            w.value  <- (17*.subset(col.work,i) + 12*(.subset(col.work,i+1)+.subset(col.work,i-1)) - 3*(.subset(col.work,i+2)+.subset(col.work,i-2))) / 35
            if(w.value<0){w.value<-0}
            col.res[i] <- w.value
          }
        }
        
      }
      
    }
    p.counts[,j]<-col.res
  }
  
  if (Round.result == T){
    p.counts <- round(p.counts)
  }
  
  return(list(Pollen=p.counts, Age=age, Age.un=data.source$Age.un, Dim.val= data.source$Dim.val))
}
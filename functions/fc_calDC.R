fc_calDC <- function (data.source, DC = "chisq")
{
  # DC = disimilarity coeficient
  #   "euc"     = Euclidan distance
  #   "euc.sd"  = standardised Euclidan distance
  #   "chord"   = chord distance
  #   "chisq    = chi-squared coeficient
  
  dat.res <- vector(mode="numeric",length = data.source$Dim.val[2]-1)
 
  
  # ----------------------------------------------
  #               EUCLIDAN DISTANCE 
  # ----------------------------------------------
  
  if (DC == "euc")
  {
    print("Euclidan distance will be used as DC")
    for (i in 1:(data.source$Dim.val[2]-1)) # for each sample (except the last)
      {
      df.work<- data.source$Pollen[c(i,i+1),] # select only 2 samples (observed + 1 after)
      
      df.work<-df.work[,colSums(df.work)>0] # get rid of "empty species"
      
      vector.work <- vector(mode="numeric", length = ncol(df.work)) # vector for result for each species
      
      for( j in 1:ncol(df.work)) # for each species
        {
        vector.work[j] <- (df.work[1,j]-df.work[2,j])**2 # calculate the diference
        }
      
      dat.res[i]<- sqrt(sum(vector.work)) # save the square root of sum of all dufereces
      
      }  
    
  }
  
  # ----------------------------------------------
  #       STANDARDISED EUCLIDAN DISTANCE 
  # ----------------------------------------------
  
  if (DC=="euc.sd")
  {
    print("Standardised Euclidan distance will be used as DC")
    
    # calculation of standard deviation for each species
    if (data.source$Dim.val[2]<1)
      stop ("too few samples for standard deviation")
    
    # dataframe for storing mean, deviance and standar deviation for each species
    df.sp.supp <- as.data.frame(matrix(ncol=2,nrow=data.source$Dim.val[1]))
    names(df.sp.supp) <- c("mean","std")
    
    for (i in 1:nrow(df.sp.supp)) # for each species
    {
      #print(paste("i",i))
      df.sp.supp$mean[i] <- mean(data.source$Pollen[,i])
      
      st.dev <- vector(mode="numeric",length =data.source$Dim.val[2])
      
      for( j in 1:data.source$Dim.val[2]) # for each sample
      {
        #print(paste("j",j))
        st.dev[j] <- (data.source$Pollen[j,i]-df.sp.supp$mean[i])**2
      }
      
      df.sp.supp$std[i]<-sqrt(sum(st.dev)/data.source$Dim.val[1])
    }
    
      # calculation of the DC
      for (i in 1:(data.source$Dim.val[2]-1)) # for each sample (except the last)
      {
        #print(paste("i",i))
        df.work<- data.source$Pollen[c(i,i+1),] # select only 2 samples (observed + 1 after)
        
        # get rid of "empty species" in data & in sp.std
        df.sp.supp.work<- df.sp.supp[colSums(df.work)>0,]
        
        df.work<-df.work[,colSums(df.work)>0] 
        
        vector.work <- vector(mode="numeric", length = ncol(df.work)) # vector for result for each species
        
        for( j in 1:ncol(df.work)) # for each species
        {
          #print(paste("j",j))
          if (df.sp.supp.work$std[j]!=0) # check if the standard deviation is not equal zero
          {
            vector.work[j] <- ((df.work[1,j]-df.work[2,j])/df.sp.supp.work$std[j])**2 # calculate the diference  
          }
          
        }
      
        dat.res[i]<- sqrt(sum(vector.work)) # save the square root of sum of all dufereces
        
      }  
      
    }
    
  
  # ----------------------------------------------
  #               CHORD DISTANCE 
  # ----------------------------------------------
  
  if (DC == "chord")
  {
    print("Chord distance will be used as DC")
    for (i in 1:(data.source$Dim.val[2]-1)) # for each sample (except the last)
    {
      df.work<- data.source$Pollen[c(i,i+1),] # select only 2 samples (observed + 1 after)
      
      df.work<-df.work[,colSums(df.work)>0] # get rid of "empty species"
      
      vector.work <- vector(mode="numeric", length = ncol(df.work)) # vector for result for each species
      
      for( j in 1:ncol(df.work)) # for each species
      {
        vector.work[j] <- (sqrt(df.work[1,j])-sqrt(df.work[2,j]))**2 # calculate the diference
      }
      
      dat.res[i]<- sqrt(sum(vector.work)) # save the square root of sum of all dufereces
      
    }  
    
  }
  
  
  # ----------------------------------------------
  #           CHI-SQUARED COEFICIENT 
  # ----------------------------------------------
  
  if (DC == "chisq")
  {
    print("Chi-squared coeficient will be used as DC")
    for (i in 1:(data.source$Dim.val[2]-1)) # for each sample (except the last)
    {
      df.work<- data.source$Pollen[c(i,i+1),] # select only 2 samples (observed + 1 after)
  
      df.work<-df.work[,colSums(df.work)>0] # get rid of "empty species"
      
      vector.work <- vector(mode="numeric", length = ncol(df.work)) # vector for result for each species
      
      for( j in 1:ncol(df.work)) # for each species
      {
        vector.work[j] <- ((df.work[1,j]-df.work[2,j])**2) / (df.work[1,j]+df.work[2,j]) # calculate the diference
      }
      
      dat.res[i]<- sqrt(sum(vector.work)) # save the square root of sum of all dufereces
      
    }  
    
  }
  
  # ----------------------------------------------
  #                   RESULT 
  # ----------------------------------------------

    return(dat.res)
  
  
}
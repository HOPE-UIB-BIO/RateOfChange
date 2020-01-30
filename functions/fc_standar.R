fc_standar <- function (data.source, fc_standar.S.value)
{
  print(paste("Data standardization started",Sys.time()))
  res.df <- data.source # prepare result df with same structure
  
  data.pol.l <- nrow(data.source) # number of samples
  
  pb<-txtProgressBar(min = 1, max=data.pol.l) # create text progres bar
  
  for(i in 1: data.pol.l)  # for each row(sample)
  {
    setTxtProgressBar(pb,i) #add progress bar unit
    temp.row <- slice(data.source,i) # select row 
    red.row <- temp.row # vector that is going to be reduced 
    res.row <- rep(0,length(temp.row)) # vector that is going to be increased
    N.row <- 1:length(temp.row) # order of vector
    
    for (j in 1:fc_standar.S.value) # Number of pollen grain we want to standardice to (150) 
    {
      N.select <- sample(N.row[red.row>0],1, prob = red.row[red.row>0]/sum(red.row)) # randomly select a species, which has more pollen than 0 
      red.row[N.select] <-  (red.row[N.select]-1) # decrease that species by 1 pollen grain in reduction vector
      res.row[N.select] <- (res.row[N.select]+1)  # increase that species by 1 pollen grain in final vector
    }
    
    res.df[i,] <- res.row # save the final vector as new values
    close(pb) # close progress bar
  }
  
  print(paste("Data standardization finished",Sys.time()))
  return (res.df) 
}
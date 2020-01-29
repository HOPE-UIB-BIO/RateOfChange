fc_standar <- function (data.source, S.value)
{
  data.pol <- data.source[[1]] # subset for only polen data
  res.df <- data.pol # prepare result df with same structure
  
  data.pol <- data.pol[,-1] # filter out the sample ID
  
  
  for(i in 1: nrow(data.pol))  # for each row
  {
    temp.row <- slice(data.pol,i) # select row 
    red.row <- temp.row # vector that is going to be reduced 
    res.row <- rep(0,length(temp.row)) # vector that is going to be increased
    N.row <- 1:length(temp.row) # order of vector
    
    for (j in 1:S.value) # Number of pollen grain we want to standardice to (150) 
    {
      N.select <- sample(N.row[red.row>0],1) # randomly select a species, which has more pollen than 0 
      red.row[N.select] <-  (red.row[N.select]-1) # decrease that species by 1 pollen grain in reduction vector
      res.row[N.select] <- (res.row[N.select]+1)  # increase that species by 1 pollen grain in final vector
    }
    
  res.df[i,-1] <- res.row # save the final vector as new values
  }
  
 return (res.df) 
}
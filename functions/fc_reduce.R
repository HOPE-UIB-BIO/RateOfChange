fc_reduce <- function (data.source.reduce,interest.treshold)
{
  data.source.reduce$Pollen <- data.source.reduce$Pollen[data.source.reduce$Age$age <= interest.treshold,]
  data.source.reduce$Age.un <- data.source.reduce$Age.un[,data.source.reduce$Age$age <= interest.treshold]
  data.source.reduce$Age <-  data.source.reduce$Age[data.source.reduce$Age$age <= interest.treshold,]
  
  
  data.source.reduce <- fc_check(data.source.reduce,proportion = F, Debug = F)
  
 return(data.source.reduce)
}
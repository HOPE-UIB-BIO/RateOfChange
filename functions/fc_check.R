fc_check <- function (data.source, proportion = F, Debug=F)
{
  # check if there is a sample that do not have a pollen data and delete it
  # & 
  #check if there are any specie without pollen record and delete them
  kill.all <- function(data.source.kill)
  {
    if(any(rowSums(data.source.kill$Pollen)==0)) 
    {
      data.source.kill$Age <- data.source.kill$Age[rowSums(data.source.kill$Pollen)>0,]
      data.source.kill$Pollen <- data.source.kill$Pollen[rowSums(data.source.kill$Pollen)>0,]
    }
    
    if(any(colSums(data.source.kill$Pollen)==0))
    {
      data.source.kill$Pollen<- data.source.kill$Pollen[colSums(data.source.kill$Pollen)>0]
    }
    
    data.source.kill$Dim.val[1]<-ncol(data.source.kill$Pollen)
    data.source.kill$Dim.val[2]<-nrow(data.source.kill$Pollen)
    data.source.kill$Dim.val[3]<-nrow(data.source.kill$Age)
    
    return(data.source.kill)
  }
  
  data.source <- kill.all(data.source) 
  
  if(Debug==T)
  {
    print("-")
    print(paste("Pollen data have",data.source$Dim.val[1],"species with pollen record and",
                data.source$Dim.val[2],"samples. Age data have",data.source$Dim.val[3],"samples"))
    print("-")
    print(paste("Age data has values of min",min(data.source$Age$age),", max",max(data.source$Age$age),",mean",
                round(mean(data.source$Age$age),2),",and median",round(median(data.source$Age$age),2)))
    print("-")
    
  }
  
  # check if all values is new age are in positive values and interpolate if necesery
  if(any(data.source$Age$newage<0))
  {
    data.source$Age$newage <- data.source$Age$newage + min(data.source$Age$newage)*(-1)
  }
  
  if (proportion == T)
  {
    if (Debug==T){print ("POllen values converted to proportions")}

    # convert the values pollen data to proportion of sum of each sample
    p.counts.row.sums <- apply(data.source$Pollen, 1, sum)
    data.source$Pollen <- as.data.frame(lapply(data.source$Pollen, function(x) x/p.counts.row.sums))
    data.source <- kill.all(data.source)
  }
  
  return(data.source)
}
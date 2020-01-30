fc_extract <-  function (data.source, standardise=T, S.value=150)
{
  # data.source = load data in format of one dataset from tibble
  # standardise = aparameter if the polen data shoudle be standardise to cetrain number of pollen grains
  # S.value = NUmber of grain to perform standardisation
  
  # result of function is list length 2
  # [1] POllen data (preferably standardise by function to S.value) saved as proportion of sum of sample
  # [2] Age data with samples ordered by age
  
  
  print(paste("Data extraction started",Sys.time()))
  
  # extract both important tables a) age data, b) pollen data
  age <- data.source$list_ages[[1]]$ages
  p.counts <- data.source$filtered.counts[[1]]
  
  # create the dimenstion of matrixes
  dim.check<-function()
  {
    assign("age.row", nrow(age), envir= parent.frame())  
    assign("p.counts.row", nrow(p.counts), envir= parent.frame())  
    assign("p.counts.col", ncol(p.counts), envir= parent.frame())  
  }
  
  #perform dimension check
  dim.check()
  
  if (age.row!=p.counts.row) # check number of rows
    stop("Pollen and Age data have different number of samples")
  
  print(paste("Pollen data have",p.counts.col,"species and",
              p.counts.row,"samples. Age data have",age.row,"samples"))
  
  # check if are sample iD saved as characters and tranform if necesery
  if(is.character(age$sample.id)==F) {age$sample.id <- as.character(age$sample.id)}  
  if(is.character(row.names(p.counts))==F) {row.names(p.counts) <- as.character(row.names(p.counts))}
  
  if (any(row.names(p.counts)!=age$sample.id)) # check if samples have same name
    stop("Samples code for pollen data and age data have different names")
  
  if(is.unsorted(age$age)==T) # check if is age of samples in order
    stop("Age data is not sorted")
  
  if(age$age[1]>age$age[age.row]){print("Age data is in decreasing format")}
  if(age$age[1]<age$age[age.row]){print("Age data is in increasing format")}
  
  if(any(rowSums(p.counts)==0)) # chek if there is a sample that do not have a pollen data and delete it
  {
    age<- age[rowSums(p.counts)>0,]
    p.counts <- p.counts[rowSums(p.counts)>0,]
  }
  
  # create a new variable for age that would be used for all latter analysys
  age$newage <-age$age
  
  # check if all values is new age are in positive values and interpolate if necesery
  if(any(age$newage<0))
  {
    age$newage <- age$newage + min(age$newage)*(-1)
  }
  
  
  if(standardise==T) # standardisation of pollen data to X(S.value) number of pollen grains
  {
    p.counts <- fc_standar(p.counts, S.value)
    if(any(rowSums(p.counts)!=S.value))
      stop("standardisation was unsuccesfull")
  }
  
  if(any(colSums(p.counts)==0)) #check if there are any specie without pollen record and delete them
  {
    p.counts<- p.counts[colSums(p.counts)>0]
  }
  
  # convert the values pollen data to proportion of sum of each sample
  p.counts.row.sums <- apply(p.counts, 1, sum)
  p.counts <- as.data.frame(lapply(p.counts, function(x) x/p.counts.row.sums))

  #perform dimension check
  dim.check()
  
  print(paste("Pollen data have",p.counts.col,"species with pollen record and",
              p.counts.row,"samples. Age data have",age.row,"samples"))
  
  print(paste("Age data has values of min",min(age$age),", max",max(age$age),",mean",
              round(mean(age$age),2),",and median",round(median(age$age),2)))
  
  
  # cretae list of 2 variables POllen & age
  dat.merge <- list(Pollen=p.counts, Age=age)
  print(paste("Data extraction completed",Sys.time()))
  return(dat.merge)
}
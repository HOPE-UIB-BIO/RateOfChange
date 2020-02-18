fc_extract <-  function (data.source.pollen.extract,
                         data.source.age.extract,
                         Debug = F)
{
  # data.source = load data in format of one dataset from tibble
  # result of function is list length 3
  # [1] POllen data 
  # [2] Age data with samples ordered by age
  # [3] NUmber of pollen species and number of samples for pollen & age data
  # [4] Dataframe with all age uncertainties
  
  if (Debug==T)
  {
    cat("",fill = T)
    cat(paste("Data extraction started",Sys.time()),fill = T)
    
  }
  
  # extract both important tables a) age data, b) pollen data
  age <- data.source.age.extract$ages
  p.counts <- data.source.pollen.extract
  age.un <- data.frame(data.source.age.extract$age_position)
  names(age.un) <- age$sample.id
  
  # create a new variableswould be used all latter analysys
  # Newgae is a value of interpolated time (time which start with 0)
  age$newage <-age$age
  
  # dim.val are values of the size of the dataset
  dim.val <- vector(mode = "integer", length = 3)
  names(dim.val) <- c("N Species","N samples pollen","N samples Age")
  
  # cretae list of 3 variables POllen, age, dim.val
  dat.merge <- list(Pollen=p.counts, Age=age, Dim.val = dim.val, Age.un = age.un)
  
  # perform check = cound number of species and samples and exclude "empty" ones
  dat.merge <- fc_check(dat.merge, proportion = F, Debug = Debug)
  
  if (dat.merge$Dim.val[3]!=dat.merge$Dim.val[2]) # check number of rows
    stop("Pollen and Age data have different number of samples")
  
  if (dat.merge$Dim.val[3]!=ncol(dat.merge$Age.un))
    stop("Age uncertainty data does not have appropriete number of samples")
  
  # check if are sample iD saved as characters and tranform if necesery
  if(is.character(dat.merge$Age$sample.id)==F) {dat.merge$Age$sample.id <- as.character(dat.merge$Age$sample.id)}  
  if(is.character(row.names(dat.merge$Pollen))==F) {row.names(dat.merge$Pollen) <- as.character(row.names(dat.merge$Pollen))}
  
  if (any(row.names(dat.merge$Pollen)!=dat.merge$Age$sample.id)) # check if samples have same name
    stop("Samples code for pollen data and age data have different names")
  
  if(is.unsorted(dat.merge$Age$age)==T) # check if is age of samples in order
    stop("Age data is not sorted")
  
  if(dat.merge$Age$age[1]>dat.merge$Age$age[dat.merge$Dim.val[3]]){
    if (Debug==T)
    {
      cat("Age data was in decreasing format, changed accordingly",fill = T)
      dat.merge$Age <- dat.merge$Age[order(dat.merge$Age$age),]
    }
    
    }
  if(dat.merge$Age$age[1]<dat.merge$Age$age[dat.merge$Dim.val[3]]){
    if (Debug==T)
    {
      cat("Age data is in increasing format",fill = T)  
    }
    
    }
  
  if (Debug==T)
  {
    cat("",fill = T)
    cat(paste("Data extraction completed",Sys.time()),fill = T)
    cat("",fill = T)
    
  }
  
  return(dat.merge)
}
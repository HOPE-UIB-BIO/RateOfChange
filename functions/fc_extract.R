fc_extract <-  function (data.source)
{
  # extract both important tables a) age data, b) pollen data
  age <- data.source$bchron_ages[[1]]$ages
  p.counts <- data.source$compiled.counts[[1]]
  
  if (nrow(age)!=nrow(p.counts)) # check number of rows
    stop("POllen and Age data have different number of rows")
  
  # check if are sample iD saved as characters and tranform if necesery
  if(is.character(age$sample.id)==F) {age$sample.id <- as.character(age$sample.id)}  
  if(is.character(p.counts$sample.id)==F) {p.counts$sample.id <- as.character(p.counts$sample.id)}
    
  # cretae list of 2 variables POllen & age
  dat.merge <- list(Pollen=p.counts, Age=age)
  return(dat.merge)
}

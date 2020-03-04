fc_extrap <- function(data.source.extrap,BIN, Debug)
{
  # function for extrapolation of new points in dataset, whie keeping all the samples
  
  # get rid of duplicate values in age#
  #AGE.DF  <- distinct(data.source.extrap$Age,newage, .keep_all = T)
  AGE.DF <- data.source.extrap$Age
  
  #ensure that there are any duplicated values in ages
  repeat {
    if(any(duplicated(AGE.DF$newage))==F){break}
    AGE.DF$newage[duplicated(AGE.DF$newage)] <- AGE.DF$newage[duplicated(AGE.DF$newage)]+1
  }
  AGE <- AGE.DF$newage

  POLEN <- data.source.extrap$Pollen
  
  # create new vector with existing samples + samples evenly distributed by BIN  
  seq.w <- seq(from= min(AGE),
               to =max(AGE),
               by = BIN) 
  seq.f <- unique(sort(c(seq.w,AGE))) 
  
  
  # create result DF
  res.DF <- data.frame(matrix(ncol=data.source.extrap$Dim.val[1],
                              nrow=length(seq.f)))
  names(res.DF) <- names(POLEN)
  
  
  # extrapolate sample ID
  AGE.res <- left_join(data.frame(newage=seq.f),AGE.DF,by="newage")
  AGE.res$sample.id <- as.numeric(AGE.res$sample.id)
  missing.ID <- c(1:nrow(AGE.res))[is.na(AGE.res$sample.id)] 
  
   
  for( j in 1:length(missing.ID))
  {
    search.number <- missing.ID[j]
    scope <- 1
      repeat
      {
        if(search.number+scope>=max(AGE.res$sample.id,na.rm = T)){break}
        if(is.na(AGE.res$sample.id[search.number+scope])!=T){break}
        scope<-scope+1
      }
      
    AGE.res$sample.id[(search.number-1):(search.number+scope)] <- seq(from=AGE.res$sample.id[(search.number-1)],
                                                                      to=AGE.res$sample.id[search.number+scope],
                                                                      length.out = scope+2)
    j<-j+(scope-1)
  }
  
  # check if every predicted value is unique
  if(length(AGE.res$sample.id) != length(unique(AGE.res$sample.id)))
  {
    # find duplicated values and add tiny value
    AGE.res$sample.id[duplicated(AGE.res$sample.id)] <- AGE.res$sample.id[duplicated(AGE.res$sample.id)]+0.0001 
  }
  
  AGE.res$sample.id <- as.character(AGE.res$sample.id)
  
  # extrapolate Pollen data
  for(i in 1:data.source.extrap$Dim.val[1]) # for each species
  {
    DF.w <- approx(x=AGE,
                 y=POLEN[,i],
                 xout = seq.f)
    res.DF[,i] <- DF.w$y
  }
  
  row.names(res.DF) <- AGE.res$sample.id
  
  list.res <- list(Pollen=res.DF, 
       Age=AGE.res,
       Dim.val=data.source.extrap$Dim.val)
  
  list.res<- fc_check(list.res, proportion = F)
  
  return(list.res)
}

fc_extrap <- function(data.source.extrap,BIN, Debug)
{
  # function for extrapolation of new points in dataset, whie keeping all the samples
  
  # get rid of duplicate values in age
  AGE.DF  <- distinct(data.source.extrap$Age,newage, .keep_all = T)
  AGE <- AGE.DF$newage

  POLEN <- data.source.extrap$Pollen
  POLEN <- POLEN[row.names(POLEN) %in% AGE.DF$sample.id,]
  
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
  extrap.samle <- approx(x=AGE.res$newage,
                         y=as.numeric(AGE.res$sample.id),
                         xout = AGE.res$newage)
  AGE.res$sample.id <- as.character(extrap.samle$y)
  
  
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

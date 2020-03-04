fc_subset_samples <- function(data.subset, BINS)
{
  data.res.age <- data.frame(matrix(ncol = ncol(data.subset$Age), nrow = nrow(BINS)))
  names(data.res.age) <- names(data.subset$Age)
  
  data.res.pollen <- data.frame(matrix(ncol=ncol(data.subset$Pollen), nrow = nrow(BINS)))
  names(data.res.pollen) <- names(data.subset$Pollen)
  
  row.names(data.res.age) <- BINS$NAME
  row.names(data.res.pollen) <- BINS$NAME
  
  BIN.size <- BINS$NAME[2]-BINS$NAME[1]
  
  for(i in 1:nrow(BINS))
  {
    selected.BIN <- BINS$NAME[i]
    
    subset.w <- data.subset$Age[data.subset$Age$newage < BINS$NAME[i]+BIN.size & 
                                  data.subset$Age$newage > BINS$NAME[i],]
    if (nrow(subset.w)>0)
    {
      subset.w$diff <- abs(subset.w$newage-selected.BIN)
      data.res.age[i,] <- subset.w[subset.w$diff==min(subset.w$diff),c(1:4)] 
      
      data.res.pollen[i,]<-  data.subset$Pollen[row.names(data.subset$Pollen) %in% data.res.age$sample.id[i],]
    }
  }

  list.res <- list(Pollen = data.res.pollen, Age=data.res.age,Dim.val=data.subset$Dim.val )
  list.res <- fc_check(list.res, proportion = F, Debug = F, Species = T, Samples = F)
  
  return(list.res)
}
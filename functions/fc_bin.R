fc_bin <- function(data.source.bin,BIN, Debug=F)
{
# BIN size determination
b.size <- BIN*2  
b.first <- 0
b.last <- ceiling(max(data.source.bin$Age$newage)/(b.size))
BIN.seq <- seq(from=b.first, to=b.last*b.size, by=b.size)

if(Debug==T)
{
  cat(paste("time window is",b.size,"and there is",b.last,"potentional bins"))
}
  
# create empty dataframes with a bin as a sample
Pollen <- data.frame(matrix(ncol = ncol(data.source.bin$Pollen),
                            nrow= length(BIN.seq)))
names(Pollen) <- names(data.source.bin$Pollen)  
row.names(Pollen) <- BIN.seq

Age <- data.frame(matrix(ncol = ncol(data.source.bin$Age),
                         nrow= length(BIN.seq)))
names(Age) <- names(data.source.bin$Age)  
Age$sample.id <- BIN.seq

  for(i in 1:length(BIN.seq)) #for each BIN
  {
    # select samples in bin length
    data.w <- dplyr::filter(data.source.bin$Age,age >= BIN.seq[i]-BIN & age <= BIN.seq[i]+BIN)
  
    if (nrow(data.w)>0) # if there are any samples selected
    {
      Age$depth[i] <- mean(data.w$depth)
      Age$age[i] <- mean(data.w$age)
      Age$newage[i] <- mean(data.w$newage)
      
      select.row <- row.names(data.source.bin$Pollen) %in% data.w$sample.id
      
      Pollen[i,]<- colSums(data.source.bin$Pollen[select.row,])  
    } 
  }

# create result list with getting rid of the empty bins 
res.list <- list(Pollen=na.omit(Pollen),
                 Age=na.omit(Age),
                 Dim.val = data.source.bin$Dim.val)

# check the data and recalculate sizes
res.list<- fc_check(res.list,proportion = F) 

return(res.list)
}

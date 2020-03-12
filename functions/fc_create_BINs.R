fc_create_BINs <-function(data.source.bin, shift.value, N.shifts)
{
  BIN <- shift.value*N.shifts
  
  BIN.last <- ceiling(max(data.source.bin$Age$newage))
  BIN.breaks <- seq(from=0, to=BIN.last, by=BIN)
  
  BIN.breaks.temp <- BIN.breaks
  BIN.breaks.fin <- vector(mode = "numeric")
  
  
  for (j in 1:N.shifts)
  {
    vector.w <- BIN.breaks.temp+shift.value*(j-1)
    BIN.breaks.fin <-c(BIN.breaks.fin,vector.w)
  }
  
  DF.sample.names <- data.frame(NAME = BIN.breaks.fin,
                                SHIFT = sort(rep(c(1:N.shifts),length(BIN.breaks))))
  
  return(DF.sample.names)
  
}


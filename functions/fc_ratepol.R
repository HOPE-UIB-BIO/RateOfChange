fc_ratepol <- function (data.source.pollen,
                        data.source.age,
                        rand = 999,
                        standardise = T, 
                        S.value = 150, 
                        sm.type = "grim", 
                        N.points = 5, 
                        range.age.max = 300, 
                        grim.N.max = 9,
                        DC = "chisq",
                        Debug = F)
{
  # data.source = data in format of one dataset from tibble
  # 
  # rand = number of randomization for estimation significance
  #
  # standardise = aparameter if the polen data shoudle be standardise to cetrain number of pollen grains
  # 
  # S.value = Number of grain to perform standardisation
  #
  # sm.type = type of smoothing applied smooting 
  #     "none"    = data will not be smoothed 
  #     "m.avg"   = moving average
  #     "grim"    = Grimm smoothing
  #     "age.w"   = age weithed 
  #     "shep"    = Shepard's 5-term filter
  #
  # N.points = Number of points for (need to be an odd number). Used for moving average, Grimm and Age-Weighted
  # grim.N.max = maximal number of samples to look in Grimm smoothing
  # range.age.max = maximal age range for both Grimm and Age-weight smoothing
  #
  # DC = disimilarity coeficient
  #   "euc"     = Euclidan distance
  #   "euc.sd"  = standardised Euclidan distance
  #   "chord"   = chord distance
  #   "chisq    = chi-squared coeficient
  #
  # Debug = show internal messages from program ? [T/F]

  
  start.time <- Sys.time()
  
  cat(paste("RATEPOL started", start.time), fill=T)
  
  
  # ----------------------------------------------
  #               DATA EXTRACTION
  # ----------------------------------------------
  
  # dataset.ID <- data.source$dataset.id
  
  # cat(paste("Data set ID",dataset.ID), fill = T)
  
  # already include data check
  data.work <- fc_extract(data.source.pollen,data.source.age, Debug=Debug) 
  
  # ----------------------------------------------
  #             RANDOMOMIZATION
  # ----------------------------------------------
  
  # detect number of cores for parallel computation
  Ncores <- detectCores()
  
  # create cluster 
  cl <- makeCluster(Ncores)
  registerDoSNOW(cl)
  
  # add all functions to the cluster
  clusterExport(cl, c("fc_standar","fc_check","fc_smooth","fc_calDC"))
  
  # create progress bar based os the number of replication
  pb <- txtProgressBar(max = rand, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  # repeat the calculation X times, whrere X is the number of randomisation
  result.tibble <- foreach(l=1:rand,.combine = rbind, .options.snow=opts) %dopar% {
    
    # ----------------------------------------------
    #             DATA STANDARFISATION
    # ----------------------------------------------
    # standardisation of pollen data to X(S.value) number of pollen grains 
    if(standardise==T) # 
    {
      data.sd <- fc_standar(data.work, S.value, Debug=Debug)
      
      if(any(rowSums(data.sd$Pollen)!=S.value))
        stop("standardisation was unsuccesfull")
    }
    
    # data check with proportioning
    data.sd.check <- fc_check(data.sd, proportion = T, Debug=Debug)
    
    # ----------------------------------------------
    #               DATA SMOOTHING
    # ----------------------------------------------
    # smooth pollen data by selected smoothing type
    data.smooth <- fc_smooth(data.sd.check, 
                             sm.type = sm.type, 
                             N.points = N.points,
                             grim.N.max = grim.N.max, 
                             range.age.max = range.age.max,
                             Debug=Debug)
    
    #data check (with proportioning ???)
    data.smooth.check <- fc_check(data.smooth, proportion = T, Debug=Debug)
    
    # ----------------------------------------------
    #               DC CALCULATION
    # ----------------------------------------------
    # calculate DC for each sample
    DC.res <- fc_calDC(data.smooth.check,DC=DC, Debug=Debug)
    
    # ----------------------------------------------
    #             AGE STANDARDISATION
    # ----------------------------------------------
    
    sample.size.work <- data.smooth.check$Dim.val[2]-1 
    
    age.diff <- vector(mode = "numeric", length = sample.size.work )
    
    for (i in 1:sample.size.work)
    {
      age.diff[i] <- abs(data.smooth.check$Age$newage[i+1]-data.smooth.check$Age$newage[i]) 
      # temporary fix for errors in age data where age difference between samples is 0
      if(age.diff[i]<1)
      {age.diff[i]<-1}
    }
    
    
    if (Debug ==T)
    {
      cat("", fill=T)
      cat(paste("The time standardisation unit (TSU) is",round(mean(age.diff),2)), fill=T)  
    }
    
  
    DC.res.s <- vector(mode = "numeric", length = sample.size.work)
    for (j in 1:sample.size.work)
    {
      DC.res.s[j] <- DC.res[j]*mean(age.diff)/age.diff[j]
    }
    
    # ----------------------------------------------
    #         RESULT OF SINGLE RAND RUN
    # ----------------------------------------------
    
      data.result <- data.frame(sample.id=data.smooth.check$Age$sample.id[1:sample.size.work], 
                                #Age=data.smooth.check$Age$newage[1:sample.size.work],
                                RoC=DC.res.s)
      
      data.result.temp <- as.data.frame(list(ID=l,DF=data.result)) 
      
      return(data.result.temp)
      # result.tibble <- rbind(result.tibble,data.result.temp)
       close(pb) # close progress bar
  }# end of the randomization
  
  # close progress bar and cluster
  close(pb)
  stopCluster(cl) 
  
  
  # ----------------------------------------------
  #             RESULTs SUMMARY
  # ----------------------------------------------
  
  # pivot all result by the randomisations
  r<- reshape2::dcast(result.tibble, DF.sample.id~ID, value.var = "DF.RoC")
  
  r.small <- r[,-1]
  
  # create new dataframe with summary of randomisation results
  r.m <-data.frame(
                   RoC.median = apply(r.small,1, FUN = function(x) median(x)),
                   RoC.se = apply(r.small,1, FUN= function(x) sd(x)/sqrt(rand)),
                   RoC.05q = apply(r.small,1, FUN= function(x) quantile(x,0.05)),
                   RoC.95q = apply(r.small,1, FUN= function(x) quantile(x,0.95))
                   )
  
  # treshold for RoC peaks is set as median of all RoC in dataset
  r.treshold <- median(r.m$RoC.median)
  
  # calculate P-value from randomisations (what number of randomisation is higher than treshold)
  #r.m$RoC.p <- apply(select(r,-c("DF.sample.id")),1,FUN=function(x) {
  #  y<- x>r.treshold
  #  y.lengt <-length(y[y==T])
  #  z <- 1-(y.lengt/rand)
  #  return(z)
  #})
  
  # mark significant peaks
  r.m$Peak <- r.m$RoC.05q>r.treshold  #r.m$RoC.p < 0.05

  # macth the samples by the sample ID
  r.m$sample.id <- r$DF.sample.id
  suppressWarnings(r.m.full <- right_join(data.work$Age,r.m, by="sample.id"))
  
  # outro
 
  end.time <- Sys.time()
  time.length <- end.time - start.time
  cat("", fill=T)
  cat(paste("RATEPOL finished", end.time,"taking",time.length, units(time.length)), fill=T)
 
 return(r.m.full)
 
}
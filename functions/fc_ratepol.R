fc_ratepol <- function (data.source,
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
  
  print (paste("RATEPOL started", start.time))
  
  
  # ----------------------------------------------
  #               DATA EXTRACTION
  # ----------------------------------------------
  dataset.ID <- data.source$dataset.id
  
  print(paste("Data set ID",dataset.ID))
  
  # already include data check
  data.work <- fc_extract(data.source, Debug=Debug) 
  
  # ----------------------------------------------
  #             RANDOMOMIZATION
  # ----------------------------------------------
  
  # save data as new dataframe
  data.sd <- data.work
  
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
      data.sd$Pollen <- fc_standar(data.work$Pollen, S.value, Debug=Debug)
      
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
      if(age.diff[i]==0)
      {age.diff[i]<-1}
    }
    
    
    if (Debug ==T)
    {
      print("-")
      print(paste("The time standardisation unit (TSU) is",round(mean(age.diff),2)))  
    }
    
    
    DC.res.s <- vector(mode = "numeric", length = sample.size.work)
    for (j in 1:sample.size.work)
    {
      DC.res.s[j] <- DC.res[j]*mean(age.diff)/age.diff[j]
    }
    
    # ----------------------------------------------
    #         RESULT OF SINGLE RAND RUN
    # ----------------------------------------------
    
      data.result <- data.frame(Age=data.smooth.check$Age$age[1:sample.size.work], RoC=DC.res.s)
      row.names(data.result) <- row.names(data.smooth.check$Pollen)[1:sample.size.work]
      
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
  r<- reshape2::dcast(result.tibble, DF.Age~ID, value.var = "DF.RoC")
  
  # create new dataframe with summary of randomisation results
  r.m <-data.frame(DF.Age=r$DF.Age,
                   RoC.mean = rowSums(select(r,-c("DF.Age"))/rand),
                   RoC.se = apply(select(r,-c("DF.Age")),1, FUN= function(x) sd(x)/sqrt(rand)),
                   RoC.05q = apply(select(r,-c("DF.Age")),1, FUN= function(x) quantile(x,0.05)),
                   RoC.95q = apply(select(r,-c("DF.Age")),1, FUN= function(x) quantile(x,0.95))
                   )
  # copy row names
  row.names(r.m) <- row.names(result.tibble[result.tibble$ID==1,])
  
  # treshold for RoC peaks is set as median of all RoC in dataset
  r.treshold <- median(r.m$RoC.mean)
  
  # calculate P-value from randomisations (what number of randomisation is higher than treshold)
  r.m$RoC.p <- apply(select(r,-c("DF.Age")),1,FUN=function(x) {
    y<- x>r.treshold
    y.lengt <-length(y[y==T])
    z <- 1-(y.lengt/rand)
    return(z)
  })
  
  # mark significant peaks
  r.m$Peak <-  r.m$RoC.p < 0.05
  
  
 end.time <- Sys.time()
 time.length <- end.time - start.time
 print("-")
 print (paste("RATEPOL finished", end.time,"taking",time.length, units(time.length)))
 
 return(list(ID=dataset.ID,Data=r.m))
 
}
fc_ratepol <- function (data.source.pollen,
                        data.source.age,
                        rand = 1000,
                        interest.treshold = F,
                        shift.value = 100,
                        N.shifts = 5,
                        standardise = T, 
                        S.value = 150, 
                        sm.type = "grim", 
                        N.points = 5, 
                        range.age.max = 500, 
                        grim.N.max = 9,
                        DC = "chisq",
                        Debug = F)
{
  # data.source = data in format of one dataset from tibble
  # 
  # rand = number of randomization for time and pollen subsample
  #
  # interest.treshold = date after wchich the data is reduced after calucaltion of RoC 
  #     but before calculation of Threshold
  #
  # BIN =
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
  data.extract <- fc_extract(data.source.pollen,data.source.age, Debug=Debug) 
  
  # ----------------------------------------------
  #               POLLEN SMOOTHING
  # ----------------------------------------------
  data.smooth <- data.extract 
  
  # smooth pollen data by selected smoothing type
  data.smooth <- fc_smooth(data.smooth, 
                           sm.type = sm.type, 
                           N.points = N.points,
                           grim.N.max = grim.N.max, 
                           range.age.max = range.age.max,
                           Debug=Debug)
  
  #data check (with proportioning ???)
  data.work <- fc_check(data.smooth, proportion = F, Debug=Debug)
  
  # ----------------------------------------------
  #               BIN creattion
  # ----------------------------------------------
  BIN.sizes <- fc_create_BINs(data.work,shift.value = shift.value, N.shifts = N.shifts)
  
  
  # ----------------------------------------------
  #             RANDOMOMIZATION
  # ----------------------------------------------
  
  # detect number of cores for parallel computation
  Ncores <- detectCores()
  
  # create cluster 
  cl <- makeCluster(Ncores)
  registerDoSNOW(cl)
  
  # add all functions to the cluster
  clusterExport(cl, c("fc_standar","fc_check","fc_smooth","fc_calDC","fc_bin","fc_extrap",
                      "filter","distinct","left_join"))
  
  # create progress bar based os the number of replication
  pb <- txtProgressBar(max = rand, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  # repeat the calculation X times, whrere X is the number of randomisation
  result.tibble <- foreach(l=1:rand,.combine = rbind, .options.snow=opts) %dopar% {
    
    # TIME SAMPLING
    # sample random time sequence from time uncern.
    data.work$Age$newage <- as.numeric(data.work$Age.un[sample(c(1:nrow(data.work$Age.un)),1),])
    
    
    # repeat for number of shifts
    for(i in 1: N.shifts)
    {
      # ----------------------------------------------
      #             DATA SUBSETTING
      # ----------------------------------------------
      SELECTED.BINS <- BIN.sizes[BIN.sizes$SHIFT==i,]
      
      data.subset <- fc_subset_samples(data.work,SELECTED.BINS)
      
      # ----------------------------------------------
      #         DATA STANDARDISATION
      # ----------------------------------------------
      data.sd <-data.subset
      
      # standardisation of pollen data to X(S.value) number of pollen grains 
      if(standardise==T) # 
      {
        data.sd <- fc_standar(data.sd, S.value, Debug=Debug)
        
        if(any(rowSums(data.sd$Pollen, na.rm = T)!=S.value & rowSums(data.sd$Pollen, na.rm = T)!=0))
          stop("standardisation was unsuccesfull")
      }
      
      # data check with proportioning
      data.sd.check <- fc_check(data.sd, proportion = T, Samples = F, Debug=Debug)
      
      
      
      
      
      
    }
    
    
    
    
    
    
    
    
    
    # ----------------------------------------------
    #               DC CALCULATION
    # ----------------------------------------------
    # calculate DC for each sample
    DC.res <- fc_calDC(data.sd.check,DC=DC, Debug=Debug)
    
    # ----------------------------------------------
    #             AGE STANDARDISATION
    # ----------------------------------------------
    
    sample.size.work <- data.sd.check$Dim.val[2]-1 
    
    age.diff <- vector(mode = "numeric", length = sample.size.work )
    
    for (i in 1:sample.size.work)
    {
      age.diff[i] <- data.sd.check$Age$newage[i+1]-data.sd.check$Age$newage[i] 
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
      DC.res.s[j] <- (DC.res[j]*mean(age.diff))/age.diff[j]
    }
    
    # ----------------------------------------------
    #         RESULT OF SINGLE RAND RUN
    # ----------------------------------------------
    
      data.result <- data.frame(sample.id=data.sd.check$Age$sample.id[1:sample.size.work],
                                Age = data.sd.check$Age$age[1:sample.size.work],
                                Newage = data.sd.check$Age$newage[1:sample.size.work],
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
  
  # create new dataframe with summary of randomisation results
  extract.res <- function (x,sel.var) 
  {
  r.m<-  select(x,c("ID","DF.sample.id",sel.var)) %>%
      pivot_wider(names_from = ID, values_from = sel.var) %>%
      mutate(
        sample.id = DF.sample.id, 
        !!sel.var := select(.,-c("DF.sample.id")) %>% apply(.,1, FUN = function(x) median(x, na.rm = T)),
        !!paste0(sel.var,".05q") := select(.,-c("DF.sample.id")) %>% apply(.,1, FUN = function(x) quantile(x,0.025, na.rm = T)),
        !!paste0(sel.var,".95q") := select(.,-c("DF.sample.id")) %>% apply(.,1, FUN = function(x) quantile(x,0.975, na.rm = T))
      ) %>%
    select(sample.id, sel.var, paste0(sel.var,".05q"), paste0(sel.var,".95q"))
    
    return(r.m)
  }

  
  # match the samples by the sample ID
  r.m.full <- right_join(extract.res(result.tibble,"DF.RoC"),
                         extract.res(result.tibble,"DF.Newage"),
                         by="sample.id")
                        
  
  # reduce results by the focus age time
  if (interest.treshold!=F)
  {
    r.m.full <- r.m.full %>%
      filter(DF.Newage<=interest.treshold)
  }
    

  # treshold for RoC peaks is set as median of all RoC in dataset
  r.treshold <- median(r.m.full$DF.RoC)
  
  # mark point which are above median
  r.m.full$soft.Peak <- r.m.full$DF.RoC > r.treshold
  
  # mark significant peaks
  r.m.full$Peak <- r.m.full$DF.RoC.05q>r.treshold  #r.m$RoC.p < 0.05

  
  # outro
 
  end.time <- Sys.time()
  time.length <- end.time - start.time
  cat("", fill=T)
  cat(paste("RATEPOL finished", end.time,"taking",time.length, units(time.length)), fill=T)
 
 return(r.m.full)
 
}
fc_ratepol <- function (data.source.pollen,
                        data.source.age,
                        sm.type = "grim", 
                        N.points = 5, 
                        range.age.max = 500, 
                        grim.N.max = 9,
                        BIN = F,
                        BIN.size = 500,
                        Shiftbin = F,
                        N.shifts = 5,
                        rand = 1000,
                        standardise = T, 
                        S.value = 150, 
                        DC = "chisq",
                        interest.treshold = F,
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
  if (BIN == F)
  {
    Shiftbin <- F
    N.shifts <- 1
  }
  
  if (BIN == T & Shiftbin == F)
  {
    N.shifts <- 1
    BIN.sizes <- fc_create_BINs(data.work,shift.value = BIN.size, N.shifts = 1)
  }
  
  if (BIN == T & Shiftbin == T)
  {
    shift.value <- BIN.size/N.shifts
    BIN.sizes <- fc_create_BINs(data.work,shift.value = shift.value, N.shifts = N.shifts)
  }
  
  
  
  # ----------------------------------------------
  #             RANDOMOMIZATION
  # ----------------------------------------------
  
  # detect number of cores for parallel computation
  Ncores <- detectCores()
  
  # create cluster 
  cl <- makeCluster(Ncores)
  registerDoSNOW(cl)
  
  # add all functions to the cluster
  clusterExport(cl, c("fc_subset_samples", "fc_standar","fc_check","fc_calDC","tibble"))
  
  # create progress bar based os the number of replication
  pb <- txtProgressBar(max = rand, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  # repeat the calculation X times, whrere X is the number of randomisation
  result.tibble <- foreach(l=1:rand,.combine = rbind, .options.snow=opts) %dopar% {
    
    # TIME SAMPLING
    # sample random time sequence from time uncern.
    data.work$Age$newage <- as.numeric(data.work$Age.un[sample(c(1:nrow(data.work$Age.un)),1),])
    
    # create result tible
    SHIFT.tibble <- tibble()
    
    # repeat for number of shifts
    for(k in 1: N.shifts)
    {
      # ----------------------------------------------
      #             DATA SUBSETTING
      # ----------------------------------------------
      data.subset <- data.work
      
      if (BIN == T)
      {
        SELECTED.BINS <- BIN.sizes[BIN.sizes$SHIFT==k,]
        
        data.subset <- fc_subset_samples(data.subset,SELECTED.BINS)  
      }
      
      # ----------------------------------------------
      #         DATA STANDARDISATION
      # ----------------------------------------------
      data.sd <-data.subset
      
      # standardisation of pollen data to X(S.value) number of pollen grains 
      if(standardise==T) # 
      {
        data.sd <- fc_standar(data.sd, S.value, Debug=Debug)
        
        if(any(rowSums(data.sd$Pollen, na.rm = T)!=S.value))
          stop("standardisation was unsuccesfull")
      }
      
      # data check with proportioning
      data.sd.check <- fc_check(data.sd, proportion = T, Samples = F, Debug=Debug)
      
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
      age.diff.names <- vector(mode = "character", length = sample.size.work )
      age.mean <- age.diff
      
      for (i in 1:sample.size.work)
      {
        age.diff[i] <- data.sd.check$Age$newage[i+1]-data.sd.check$Age$newage[i] 
        # temporary fix for errors in age data where age difference between samples is 0
        if(age.diff[i]<1)
        {age.diff[i]<-1}
        
        #calculate the average position 
        age.mean[i] <- mean(c(data.sd.check$Age$newage[i+1],data.sd.check$Age$newage[i]))
        
        # create vector with bin names
        age.diff.names[i] <- paste(row.names(data.sd.check$Age)[i],"-",row.names(data.sd.check$Age)[i+1])
      }
      
      if (Debug ==T)
      {
        cat("", fill=T)
        cat(paste("The time standardisation unit (TSU) is",round(mean(age.diff),2)), fill=T)  
      }
      
      #  calculate DC standardise by time
      DC.res.s <- vector(mode = "numeric", length = sample.size.work)
      DC.res.s <- (DC.res*mean(age.diff))/age.diff
      
      # ----------------------------------------------
      #             Result of single SHIFT
      # ----------------------------------------------
      SHIFT.tibble <- rbind(SHIFT.tibble,
                            data.frame(BIN= age.diff.names,
                                     DC = DC.res,
                                     Age.Pos = age.mean,
                                     Age.Diff = age.diff,
                                     RoC = DC.res.s,
                                     SHIFT = rep(k,sample.size.work))
                            )
    }
    
    # ----------------------------------------------
    #         RESULT OF SINGLE RAND RUN
    # ----------------------------------------------
      
      data.result.temp <- as.data.frame(list(ID=l,RUN=SHIFT.tibble)) 
      
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
  r.m<-  select(x,c("ID","RUN.BIN",sel.var)) %>%
      pivot_wider(names_from = "ID", values_from = sel.var) %>%
    mutate(BIN = sub(" -.*","",RUN.BIN) %>% as.numeric()) %>%
    arrange(.,BIN) %>%
    select(-c(BIN))
    
  N.notNA<-  r.m  %>% select(-c(RUN.BIN)) %>%
    apply(.,1, FUN = function(x){
      y<- is.na(x) %>% table ()
      return(y[1]) 
    }) 
  r.m.sel <- r.m[N.notNA>0.1*rand,]

  r.m.sum <- r.m.sel %>% 
    mutate(
        sample.id = RUN.BIN, 
        !!sel.var := select(.,-c("RUN.BIN")) %>% apply(.,1, FUN = function(x) median(x, na.rm = T)),
        !!paste0(sel.var,".05q") := select(.,-c("RUN.BIN")) %>% apply(.,1, FUN = function(x) quantile(x,0.025, na.rm = T)),
        !!paste0(sel.var,".95q") := select(.,-c("RUN.BIN")) %>% apply(.,1, FUN = function(x) quantile(x,0.975, na.rm = T))
      ) %>%
    select(sample.id, sel.var, paste0(sel.var,".05q"), paste0(sel.var,".95q"))
    
    return(r.m.sum)
  }

  # match the samples by the sample ID
  r.m.full <- right_join(extract.res(result.tibble,"RUN.RoC"),
                         extract.res(result.tibble,"RUN.Age.Pos"),
                         by="sample.id")
                        
  
  # reduce results by the focus age time
  if (interest.treshold!=F)
  {
    r.m.full <- r.m.full %>%
      filter(RUN.Age.Pos<=interest.treshold)
  }
    

  # treshold for RoC peaks is set as median of all RoC in dataset
  r.treshold <- median(r.m.full$RUN.RoC)
  
  # mark point which are above median
  r.m.full$soft.Peak <- r.m.full$RUN.RoC > r.treshold
  
  # mark significant peaks
  r.m.full$Peak <- r.m.full$RUN.RoC.05q>r.treshold  #r.m$RoC.p < 0.05

  # mark points that are abowe the GAM model (exactly 2*SD higher than GAM prediction)
  pred.gam <-  predict.gam(gam(RUN.RoC~s(RUN.Age.Pos), data = r.m.full))
  pred.gam.diff <- r.m.full$RUN.RoC - pred.gam
  r.m.full$Peak.gam <- pred.gam.diff > 2*sd(pred.gam.diff)
  
  # outro
 
  end.time <- Sys.time()
  time.length <- end.time - start.time
  cat("", fill=T)
  cat(paste("RATEPOL finished", end.time,"taking",time.length, units(time.length)), fill=T)
 
 return(r.m.full)
 
}
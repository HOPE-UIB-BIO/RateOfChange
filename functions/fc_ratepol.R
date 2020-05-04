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
  # smoothing of the pollen data
  # sm.type = type of smoothing applied for the each of the pollen type 
  #     "none"    = data will not be smoothed 
  #     "m.avg"   = moving average
  #     "grim"    = Grimm smoothing
  #     "age.w"   = age weithed 
  #     "shep"    = Shepard's 5-term filter
  # N.points = Number of points for (need to be an odd number). Used for moving average, Grimm and Age-Weighted
  # grim.N.max = maximal number of samples to look in Grimm smoothing
  # range.age.max = maximal age range for both Grimm and Age-weight smoothing
  #
  # BIN [T/F] =  subsert the data into bins of various size. In each bin, one sample is selected (closest to the starting point of the bin),
  #                 rest is discardet. Bins without any samples are also dicarded.
  # BIN.size = size of the bin (in years)
  #
  # Shiftbin [T/F] = setting allowing subseting the BIN into number of Shifts. Anaylsis is rerun with same BIN size
  #                     but the begining of each BIN is moved by BIN.size/N.shifts. 
  # N.shifts = value determining the number of shifts in BIN
  #
  # rand = number of randomization. Age sequence is randomly sampled from age-depth model uncertainties at the begining of each run. 
  #
  # standardise [T/F] = standardise each sample/BIN to cetrain number of pollen grains (random resampling without repetition)
  # S.value = Number of grain to perform standardisation
  #
  # DC = disimilarity coeficient. Type of calculation of differences between samples/BINs
  #   "euc"     = Euclidan distance
  #   "euc.sd"  = standardised Euclidan distance
  #   "chord"   = chord distance
  #   "chisq    = chi-squared coeficient
  #
  # Comment 1:
  # Note that DC is calculated between each subsequent samples/BINs and then standardise by the time difference between samples/BINs. Age of each RoC values
  #   is determined as mean between values from samples/BINs the RoC was calculated.
  #
  # Comment 2:
  # Due to randomisation, there is a chanche than there will be different combination of samples in each run (some samples/BINs might "fallout" due to 
  #     random time sampling and pollen subsapling). Each Sample-Sample (or BIN-BIN) combination will be saved and the point (combination) is included only if
  #     there was at least 0.1*rand number (10%) of results for that point (combination). RoC and time for each point (combination) is calculated as
  #     median of all results from randomisations
  #
  # interest.treshold [T/F] = age after which is the data reduced after calucaltion of RoC 
  # 
  # Debug [T/F] = show messages from internal processes
  # 
  # 
  # Significance Comment:
  # Each RoC point is validated in three ways:
  #     1) soft Peaks = Treshold value is set for whole dataset (after subseting for interest.treshold) as median of all RoC values.
  #                     Point is consider significat if it is higher than Treshold
  #     2) Peaks      = Treshold is set same as in 1). Peak is  consider significant if 95% quantile (gain from randomisation) is higher than treshold
  #     3) GAM        = Gam model is fitted with RoC and Age. Differences between GAM and each point is calculated. SD is calculated from all the differences
  #                     Peak is considered significat if it is 1.5 SD higher than GAM. 
  
  # CODE
  start.time <- Sys.time()
  
  cat(paste("RATEPOL started", start.time), fill=T)
  
  
  # ----------------------------------------------
  #               DATA EXTRACTION
  # ----------------------------------------------
  
  # extract data into working format
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
  
  #check data and reduce data dimentions  
  data.work <- fc_check(data.smooth, proportion = F, Debug=Debug)
  
  # ----------------------------------------------
  #               BIN creattion
  # ----------------------------------------------

  # series of checks for BIN setting
  if (BIN == F) # do not run Shifts if BIN is FALSE
  {
    Shiftbin <- F
    N.shifts <- 1
  }
  
  # create a list of BINs 
  if (BIN == T & Shiftbin == F)
  {
    N.shifts <- 1
    BIN.sizes <- fc_create_BINs(data.work,shift.value = BIN.size, N.shifts = 1)
  }
  
  # Calculate the shift value from BIN and N.shifts and then calculate ALL BINs (including shift)
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
  clusterExport(cl, c("fc_subset_samples","fc_standar","fc_check","fc_calDC","tibble"))
  
  # create progress bar based os the number of replication
  pb <- txtProgressBar(max = rand, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  # repeat the calculation X times, whrere X is the number of randomisation
  result.tibble <- foreach(l=1:rand,.combine = rbind, .options.snow=opts) %dopar% {
    
    # TIME SAMPLING
    # sample random time sequence from time uncern.
    data.work$Age$newage <- as.numeric(data.work$Age.un[sample(c(1:max(1,nrow(data.work$Age.un))),1),])
    
    # create result tible
    SHIFT.tibble <- tibble()
    
    # repeat for number of shifts
    for(k in 1: N.shifts)
    {
      # ----------------------------------------------
      #             DATA SUBSETTING
      # ----------------------------------------------
      data.subset <- data.work
     
      # select one sample for each bin based on the age of the samples. Sample is chones if it is the closes one to the upper end of the BIN 
      if (BIN == T)
      {
        # select BIN for this shift
        SELECTED.BINS <- BIN.sizes[BIN.sizes$SHIFT==k,]
        
        
        #subset data
        data.subset <- fc_subset_samples(data.subset,SELECTED.BINS)  
      }
      
      # ----------------------------------------------
      #         DATA STANDARDISATION
      # ----------------------------------------------
      data.sd <-data.subset
      
      # standardisation of pollen data to X(S.value) number of pollen grains 
      if(standardise==T) # 
      {
        # check if all samples has S.value of pollen grains
        data.sd$Age <- data.sd$Age[rowSums(data.sd$Pollen, na.rm = T)>=S.value,]
        data.sd$Age.un <- data.sd$Age.un[,rowSums(data.sd$Pollen, na.rm = T)>=S.value]
        data.sd$Pollen <- data.sd$Pollen[rowSums(data.sd$Pollen, na.rm = T)>=S.value,]
        data.sd<- fc_check(data.sd, proportion = F, Samples = T, Debug=Debug)
        
        
        # standardisation
        data.sd <- fc_standar(data.sd, S.value, Debug=Debug)
        
        if(any(rowSums(data.sd$Pollen, na.rm = T)!=S.value))
          stop("standardisation was unsuccesfull")
      }
      
      # data check with proportioning
      data.sd.check <- fc_check(data.sd, proportion = T, Samples = F, Debug=Debug)
      
      # ----------------------------------------------
      #               DC CALCULATION
      # ----------------------------------------------
      
      # calculate DC between each subsequent samples/BINs
      DC.res <- fc_calDC(data.sd.check,DC=DC, Debug=Debug)
      
      # ----------------------------------------------
      #             AGE STANDARDISATION
      # ----------------------------------------------
      
      # create empty vector with size = numeber of samples-1
      sample.size.work <- data.sd.check$Dim.val[2]-1 
      
      # create empty vectors for age difference calcualtion
      age.diff <- vector(mode = "numeric", length = sample.size.work )
      age.diff.names <- vector(mode = "character", length = sample.size.work )
      age.mean <- age.diff
      
      for (i in 1:sample.size.work) # for each RoC
      {
        # calcualte the age difference between subsequesnt samples
        age.diff[i] <- data.sd.check$Age$newage[i+1]-data.sd.check$Age$newage[i] 
        
        # Set age difference as 1, if age difference between samples is smaller than 1
        if(age.diff[i]<1)
        {age.diff[i]<-1}
        
        #calculate the average position of RoC
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
      
      # add the results from this shift into the result tibble
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
      
      # save result from single randomisation into data.frame with number of randomisation as ID.
      data.result.temp <- as.data.frame(list(ID=l,RUN=SHIFT.tibble)) 
      
      return(data.result.temp)
       close(pb) # close progress bar
  }# end of the randomization
  
  # close progress bar and cluster
  close(pb)
  stopCluster(cl) 
  
  
  # ----------------------------------------------
  #             RESULTs SUMMARY
  # ----------------------------------------------
  
  # create new dataframe with summary of randomisation results
  
  #function for extracting results 
  extract.res <- function (x,sel.var) 
  {
  
  # cretae pivot table randomisation ID by BIN code and order it by BIN codes  
  r.m<-  select(x,c("ID","RUN.BIN",sel.var)) %>% 
      pivot_wider(names_from = "ID", values_from = sel.var) %>%
    mutate(BIN = sub(" -.*","",RUN.BIN) %>% as.numeric()) %>%
    arrange(.,BIN) %>%
    select(-c(BIN))
  
  # calculate number of NOT NA values for each BIN code
  N.notNA<-  r.m  %>% select(-c(RUN.BIN)) %>%
    apply(.,1, FUN = function(x){
      y<- is.na(x) %>% table ()
      return(y[1]) 
    }) 
  
  # include only bin codes if they have ast least 10% of number of randomisation 
  r.m.sel <- r.m[N.notNA>0.1*rand,]

  # calculate median and 95% quantile for selected variable for each sample (BIN)
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

  # extract results and match them by BIN
  r.m.full <- right_join(extract.res(result.tibble,"RUN.RoC"),
                         extract.res(result.tibble,"RUN.Age.Pos"),
                         by="sample.id")
                        
  
  # reduce results by the focus age time
  if (interest.treshold!=F)
  {
    r.m.full <- r.m.full %>%
      filter(RUN.Age.Pos<=interest.treshold)
  }
  
  # sort samples by age
  r.m.full <- r.m.full %>%
    arrange(RUN.Age.Pos)
    
  # Median peak treshold
  # treshold for RoC peaks is set as median of all RoC in dataset
  r.treshold <- median(r.m.full$RUN.RoC)
  # mark point which are above the treshold as Peak.treshold
  r.m.full$Peak.treshold <- r.m.full$RUN.RoC > r.treshold
  # mark peaks which have 95% quantile above the treshold as Peak.treshold.95
  r.m.full$Peak.treshold.95 <- r.m.full$RUN.RoC.05q>r.treshold  

  # GAM
  # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
  pred.gam <-  predict.gam(gam(RUN.RoC~s(RUN.Age.Pos), data = r.m.full))
  pred.gam.diff <- r.m.full$RUN.RoC - pred.gam
  r.m.full$Peak.gam <- pred.gam.diff > 1.5*sd(pred.gam.diff)
  
  # SNI  
  # set moving window of 5 times higher than average distance between samples
  mean.age.window <- 5 * mean( diff(r.m.full$RUN.Age.Pos) )
  # calculate SNI (singal to noise ratio)
  SNI.calc <- CharSNI(data.frame(r.m.full$RUN.Age.Pos, r.m.full$RUN.RoC, pred.gam),mean.age.window)
  # mark points with SNI higher than 3
  r.m.full$Peak.SNI <- SNI.calc$SNI > 3
  
  
  
  # outro
 
  end.time <- Sys.time()
  time.length <- end.time - start.time
  cat("", fill=T)
  cat(paste("RATEPOL finished", end.time,"taking",time.length, units(time.length)), fill=T)
 
 return(r.m.full)
 
}
#end of CODE
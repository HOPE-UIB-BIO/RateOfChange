fc_R_ratepol <- function (data.source.pollen,
                          data.source.age,
                          sm.type = "grim", 
                          N.points = 5, 
                          range.age.max = 500, 
                          grim.N.max = 9,
                          Working.Unit = "MW",
                          BIN.size = 500,
                          N.shifts = 5,
                          rand = 1000,
                          standardise = T, 
                          S.value = 150, 
                          DC = "Chord",
                          interest.treshold = F,
                          Peak = "GAM",
                          Debug = F)
{
  # data.source.pollen = pollen data with species as collumns and levels as rows, level ID as row names  
  # data.source.age = list of 2: 
  #                       $ages = datframe
  #                                   $sample.id = unique ID of each level
  #                                   #age = age of level
  #                       $age_position = matrix with number of collumns as number of levels. Each column is one level
  #                                         each row is one age sequence from bchron
  #                                         
  # smoothing of the pollen data
  # sm.type = type of smoothing applied for the each of the pollen type 
  #     "none"    = None: Pollen data is not smoothed
  #     "m.avg"   = Moving average: Each focus value is calculated as average over N(N.points) number of levels 
  #                   (preferably ½ N before and ½ after focus level, those values are adjusted in the beginning 
  #                   and end of the core. N must be odd number). Note that each calculation is done from scratch 
  #                   and results are saved separately in order to avoid cumulative rounding errors. 
  #     "grim"    = 	Grimm’s smoothing: Similar to moving average but N is not fixed. For each level, N is selected 
  #                   as odd number between N_a (N.points) and N_b (grim.N.max), while the maintaining the maximum age 
  #                   difference from the selected levels as range.age.max. 
  #     "age.w"   = 	Age weighted average:  Similar to moving average but average is weighted by the age difference 
  #                   from the focus level and multiplied by 1/range.age.max. To avoid up-weighting levels, 
  #                   if range.age.max/AGEDIFF exceeds 1, it is saved as 1. This means that levels closer than 
  #                   range.age.max to the target age are given full weighting, but those farther away are downweighed 
  #                   by an amount increasing with age difference.
  #     "shep"    =   Shepard's 5-term filter: Smoothing over 5 points following equation: 
  #                   V_NEW=(17*V + 12*(V_((+1) )  + V_((-1) ) )-3*(V_((+2) )  + V_((-2) )))/35 ,
  #                   where V is focal level value. All values that result smaller than zero are saved as zero
  # N.points = Number of points for (need to be an odd number). Used for moving average, Grimm and Age-Weighted
  # grim.N.max = maximal number of samples to look in Grimm smoothing
  # range.age.max = maximal age range for both Grimm and Age-weight smoothing
  #
  # Working.Unit = selection of units that the DC will be calculated between
  #     "levels"  = DC is calculated between all subssequent levels
  #     "BINs"    = Selective Binning: BINs (age brackets of various size (BIN.size)) are created and  one level is selected as a representation of each BIN.
  #                 BINS serve as Working Units. BINs without any pollen data are excluded
  #     "MU"      = Selective Binning with Moving Window: This method follows a simple sequence: BINs are created, levels are selected, and RoC between BINs 
  #                 is calculated, similar to Selective Binning. However, the brackets of BINs (window) are then moved forward by selected amount
  #                 of time (Z), levels are selected again, and RoC calculated for a new set of WU. This is repeated N.shifts times 
  #                 (where Z = BIN.size/N.shifts) while keeping all the results.
  # BIN.size = size of the BIN (in years)
  # N.shifts = value determining the number of shifts of window
  #
  # rand = number of randomization. Age sequence is randomly sampled from age-depth model uncertainties 
  #         at the begining of each run. 
  #
  # standardise [T/F] = standardise each Working Unit to cetrain number of pollen grains (random resampling without repetition)
  # S.value = Number of grain to perform standardisation to
  #
  # DC = disimilarity coeficient. Type of calculation of differences between samples/BINs
  #   "euc"     = 	Euclidean distance: √(∑_(i=1)^n〖(A_i-B_i)〗^2 ), where A_i and B_i are values for Working Unit A and B given for species i.
  #   "euc.sd"  = 	Standardised Euclidean distance: √(∑_(i=1)^n〖((A_i-B_i)/〖SD〗_i )〗^2 ) , 
  #                 where A_i and B_i are values for Working Unit A and B given for species i,and 〖SD〗_i is a standard deviation for species i,
  #                 calculated from whole sequence.
  #   "chord"   = 	Chord distance:√(∑_(i=1)^n〖(√(A_i )-√(B_i ))〗^2 ), where A_i and B_i are values for Working Unit A and B given for species i.
  #   "chisq    = 	Chi-squared coefficient:√(∑_(i=1)^n〖(A_i-B_i)〗^2/((A_i+B_i))), where A_i and B_i are values for Working Unit A and B 
  #                 given for species i
  # 
  # Comment 1:
  # Note that DC is calculated between each subsequent Working Units and then standardise by the time difference between them. Age of each RoC values
  #   is determined as mean between values of Wokring Units.
  #
  # Comment 2:
  # Due to randomisation, there is a chanche than there will be different Working Unit combination  in each run (some Working Units might "fallout" due to 
  #     random time sampling and pollen subsampling). The calculation between two subsequent WU (i.e. one Working Unit combination) results in RoC 
  #     score and time position (calculated as mean time position of WUs). However, due to random selecting of age sequence, each of the Working Unit
  #     combination will result in multiple RoC values and time positions. R-Ratepol assign the time position of each Working Unit combination as median
  #     time positions from all randomisations. Final RoC values are calculated as median of scores from all randomisations. In addition, due to 
  #     excluding empty Working Units, there is a chance that some Working Unit combination will be present only in some randomisation. 
  #     Therefore, R-Ratepol only include Working Unit combinations that are present in at least 10% of all randomisations.  
  #
  # interest.treshold [T/F] = age after which is the data reduced after calucaltion of RoC 
  # 
  # Peak = method of peak detection:
  #         1) Threshold = Treshold value is set for whole dataset (after subseting for interest.treshold) as median of all RoC values.
  #                     Peak is  consider significant if 95% quantile (gain from randomisation) is higher than treshold
  #         2) GAM = Gam model is fitted with RoC and Age. Differences between GAM and each point is calculated. SD is calculated from all the differences
  #                     Peak is considered significat if it is 1.5 SD higher than GAM
  #         3) SNI = Signal-to-Noise Index, following adapted the SNI from Kelly et al. (2011) 
  #                     written to detect changes in charcoal records. We SNI I calculated for the whole RoC sequence 
  #                     and point is consider significant if has SNI value higher than 3 (following suggesting 
  #                     from Kelly et al. (2011)).   
  #
  # Debug [T/F] = show messages from internal processes
   
  
  # CODE
  start.time <- Sys.time()
  
  cat(paste("RATEPOL started", start.time), fill=T)
  
  
  # ----------------------------------------------
  #               DATA EXTRACTION
  # ----------------------------------------------
  
  # extract data into working format
  # already include data check
  data.extract <- fc_extract_data(data.source.pollen,data.source.age, Debug=Debug) 
  
  # ----------------------------------------------
  #               POLLEN SMOOTHING
  # ----------------------------------------------
  data.smooth <- data.extract 
  
  # smooth pollen data by selected smoothing type
  data.smooth <- fc_smooth_pollen_data(data.smooth, 
                                        sm.type = sm.type, 
                                        N.points = N.points,
                                        grim.N.max = grim.N.max, 
                                        range.age.max = range.age.max,
                                       Round.result = standardise,
                                        Debug=Debug)
  
  #check data and reduce data dimentions  
  data.work <- fc_check_data(data.smooth, proportion = F, Debug=Debug)
  
  # ----------------------------------------------
  #               Working Unit selection 
  # ----------------------------------------------

  if(Working.Unit == "levels"){
    N.shifts <- 1
  }
  
  if(Working.Unit == "BINs"){
    N.shifts <- 1
    BIN.sizes <- fc_create_BINs(data.work,shift.value = BIN.size, N.shifts = 1)
  }
  
  if(Working.Unit == "MW"){
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
  clusterExport(cl, c("fc_subset_samples","fc_standardise_pollen_data","fc_check_data","fc_calculate_DC","tibble"))
  
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
      if (Working.Unit != "levels")
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
        
        # adjust the Svalue by the minimal Pollen or to a minimal of presected values
        S.value <-  min(c( rowSums(data.subset$Pollen),S.value))
        
        # check if all samples has S.value of pollen grains
        data.sd$Age <- data.sd$Age[rowSums(data.sd$Pollen, na.rm = T)>=S.value,]
        data.sd$Age.un <- data.sd$Age.un[,rowSums(data.sd$Pollen, na.rm = T)>=S.value]
        data.sd$Pollen <- data.sd$Pollen[rowSums(data.sd$Pollen, na.rm = T)>=S.value,]
        data.sd<- fc_check_data(data.sd, proportion = F, Samples = T, Debug=Debug)
        
        # standardisation
        data.sd <- fc_standardise_pollen_data(data.sd, S.value, Debug=Debug)
        
        if(any(rowSums(data.sd$Pollen, na.rm = T)!=S.value))
          stop("standardisation was unsuccesfull")
      }
      
      # data check with proportioning
      data.sd.check <- fc_check_data(data.sd, proportion = T, Samples = F, Debug=Debug)
      
      # ----------------------------------------------
      #               DC CALCULATION
      # ----------------------------------------------
      
      # calculate DC between each subsequent samples/BINs
      DC.res <- fc_calculate_DC(data.sd.check,DC=DC, Debug=Debug)
      
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
  r.m<-  x %>%
    mutate(BIN_shift = paste0(RUN.SHIFT,"_",RUN.BIN)) %>%
    dplyr::select(c(ID,RUN.BIN,RUN.SHIFT,sel.var)) %>% 
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
  r.m.sel <- r.m %>%
    filter(N.notNA>0.1*rand)

  # calculate median and 95% quantile for selected variable for each sample (BIN)
  r.m.sum <- r.m.sel %>% 
    mutate(
        sample.id = RUN.BIN, 
        shift = RUN.SHIFT,
        !!sel.var := select(.,-c(RUN.BIN,RUN.SHIFT)) %>% apply(.,1, FUN = function(x) median(x, na.rm = T)),
        !!paste0(sel.var,".05q") := select(.,-c(RUN.BIN,RUN.SHIFT)) %>% apply(.,1, FUN = function(x) quantile(x,0.025, na.rm = T)),
        !!paste0(sel.var,".95q") := select(.,-c(RUN.BIN,RUN.SHIFT)) %>% apply(.,1, FUN = function(x) quantile(x,0.975, na.rm = T))
      ) %>%
    select(sample.id, shift, sel.var, paste0(sel.var,".05q"), paste0(sel.var,".95q"))
    
    return(r.m.sum)
  }

  # extract results and match them by BIN
  r.m.full <- right_join(extract.res(result.tibble,"RUN.RoC"),
                         extract.res(result.tibble,"RUN.Age.Pos") %>%
                           select(-shift),
                         by="sample.id")
                        
  
  # reduce results by the focus age time
  if (interest.treshold!=F)
  {
    r.m.full <- r.m.full %>%
      filter(RUN.Age.Pos<=interest.treshold)
  }
  
  # sort samples by age and add smoothing to avoid "waveing"
  r.m.full <- r.m.full %>%
    arrange(RUN.Age.Pos) %>%
    mutate (
      RUN.RoC.sm =  lowess(RUN.Age.Pos,RUN.RoC,f=.1,iter=100)$y,
      RUN.RoC.95q.sm = lowess(RUN.Age.Pos,RUN.RoC.95q,f=.1,iter=100)$y,
      RUN.RoC.05q.sm = lowess(RUN.Age.Pos,RUN.RoC.05q,f=.1,iter=100)$y
    ) %>%
    mutate(RUN.RoC.sm = ifelse(RUN.RoC.sm<0,0,RUN.RoC.sm),
           RUN.RoC.05q.sm = ifelse(RUN.RoC.05q.sm<0,0,RUN.RoC.05q.sm))
  
  
  # ----------------------------------------------
  #             PEAK DETECTION
  # ----------------------------------------------
  
  # Median peak treshold
  if(Peak == "Threshold"){
    # treshold for RoC peaks is set as median of all RoC in dataset
    r.treshold <- median(r.m.full$RUN.RoC.sm)
    # mark peaks which have 95% quantile above the treshold asPeak.treshold
    r.m.full<- r.m.full %>%
      mutate(Peak = RUN.RoC.05q.sm > r.treshold)
  }
  
  # GAM  
  if(Peak == "GAM"){
    # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
    r.m.full<- r.m.full %>%
      mutate(
        pred.gam = predict.gam(gam(RUN.RoC.sm~s(RUN.Age.Pos,k=3), data = ., family = "Gamma",
                                   method = "REML"), type="response"),
        pred.gam.diff = RUN.RoC.sm - pred.gam,
        Peak = (pred.gam.diff) > 1.5*sd(pred.gam.diff)
        )
  }
  
  # SNI  
  if (Peak == "SNI"){
    # set moving window of 5 times higher than average distance between samples
    mean.age.window <- 5 * mean( diff(r.m.full$RUN.Age.Pos) )
    # create GAM 
    pred.gam <-  predict.gam(gam(RUN.RoC.sm~s(RUN.Age.Pos,k=3), data = r.m.full, family = "Gamma",
                                 method = "REML"), type="response")
    # calculate SNI (singal to noise ratio)
    SNI.calc <- CharSNI(data.frame(r.m.full$RUN.Age.Pos, r.m.full$RUN.RoC.sm, pred.gam),mean.age.window)
    # mark points with SNI higher than 3
    r.m.full <- r.m.full %>%
      mutate( Peak = SNI.calc$SNI > 3 & r.m.full$RUN.RoC.sm > pred.gam )
  }
  
  # outro
 
  r.m.full.fin <- r.m.full %>%
    rename(Working_Unit = sample.id,
           ROC = RUN.RoC.sm,
           ROC.up = RUN.RoC.95q.sm,
           ROC.dw =RUN.RoC.05q.sm,
           AGE = RUN.Age.Pos,
           PEAK = Peak) %>%
    dplyr::select(Working_Unit,ROC,ROC.up,ROC.dw,AGE,PEAK)
  
  
  
  end.time <- Sys.time()
  time.length <- end.time - start.time
  cat("", fill=T)
  cat(paste("RATEPOL finished", end.time,"taking",time.length, units(time.length)), fill=T)
 
  
 return(r.m.full.fin)
 
}
#end of CODE
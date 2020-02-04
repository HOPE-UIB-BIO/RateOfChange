fc_ratepol <- function (data.source, 
                        standardise = T, 
                        S.value = 150, 
                        sm.type = "grim", 
                        N.points = 5, 
                        range.age.max = 300, 
                        grim.N.max = 9,
                        DC = "chisq",
                        result = "tibble")
{
  # data.source = data in format of one dataset from tibble
  # standardise = aparameter if the polen data shoudle be standardise to cetrain number of pollen grains
  # S.value = NUmber of grain to perform standardisation
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
  # result = type of requested result
  #   "small" = result list of 2: [1] DF with age and Rate of change; 
  #                               [2] plot of Rate of change in time
  #   "full"  = result list of 2: [1] DF with all pollen samples, age data, Rate of Change data, 
  #                               [2] plot of species pollen changes and Rate of Change in time 
  
  start.time <- Sys.time()
  
  print("-")
  print (paste("RATEPOL started", start.time))
  
  
  # ----------------------------------------------
  #               DATA EXTRACTION
  # ----------------------------------------------
  dataset.ID <- data.source$dataset.id
  print(paste("Data set ID",dataset.ID))
  print("-")
  
  # already include data check
  data.work <- fc_extract(data.source) 
  
  # ----------------------------------------------
  #             DATA STANDARFISATION
  # ----------------------------------------------
  # standardisation of pollen data to X(S.value) number of pollen grains 
  if(standardise==T) # 
  {
    data.work$Pollen <- fc_standar(data.work$Pollen, S.value)
    
    if(any(rowSums(data.work$Pollen)!=S.value))
      stop("standardisation was unsuccesfull")
  }
  
  # data check with proportioning
  data.work <- fc_check(data.work, proportion = T)
  
  # ----------------------------------------------
  #               DATA SMOOTHING
  # ----------------------------------------------
  data.smooth <- fc_smooth(data.work, 
                           sm.type = sm.type, 
                           N.points = N.points,
                           grim.N.max = grim.N.max, 
                           range.age.max = range.age.max)
  
  #data check (with proportioning ???)
  data.smooth <- fc_check(data.smooth, proportion = T)
  
  # ----------------------------------------------
  #               DC CALCULATION
  # ----------------------------------------------
  # calculate DC for each sample
  DC.res <- fc_calDC(data.smooth,DC=DC)
  
  # ----------------------------------------------
  #             AGE STANDARDISATION
  # ----------------------------------------------
  
  sample.size.work <- data.smooth$Dim.val[2]-1 
  
  age.diff <- vector(mode = "numeric", length = sample.size.work )
  for (i in 1:sample.size.work)
  {
    age.diff[i] <- abs(data.smooth$Age$newage[i+1]-data.smooth$Age$newage[i])  
  }
  
  print("-")
  print(paste("The time standardisation unit (TSU) is",round(mean(age.diff),2)))
  
  DC.res.s <- vector(mode = "numeric", length = sample.size.work)
  for (j in 1:sample.size.work)
  {
    DC.res.s[j] <- DC.res[j]*mean(age.diff)/age.diff[j]
  }
  
  # ----------------------------------------------
  #               RESULT SUMMARY
  # ----------------------------------------------
  
  
  if (result == "full")
  {
    data.plot <- data.frame(data.smooth$Pollen[1:sample.size.work,],RoC=DC.res.s)
    
    data.plot.melt <- reshape2::melt(data.plot)
    data.plot.melt$age <- c(rep(data.smooth$Age$age[1:sample.size.work],data.smooth$Dim.val[1]+1)) 
    
    p1<-ggplot(data = data.plot.melt, 
               aes(y=value, 
                   x= age,
                   color=variable))+
      theme_classic()+
      scale_x_continuous(trans = "reverse")+
      coord_flip(ylim = c(0,1))+
      #geom_point(alpha=1/5)+
      geom_line()
    
    data.result <- data.frame(data.smooth$Pollen,
                              data.smooth$Age,
                              Age.diff=c(age.diff,NA),
                              DC=c(DC.res,NA),
                              Roc=c(DC.res.s,NA))
    
    row.names(data.result) <- row.names(data.smooth$Pollen)
    r<-list(ID=dataset.ID, data=data.result, plot=p1)
  }
  
  if (result == "small")
  {
    data.result <- data.frame(Age=data.smooth$Age$age[1:sample.size.work], RoC=DC.res.s)
    row.names(data.result) <- row.names(data.smooth$Pollen)[1:sample.size.work]
    
    p1<- ggplot(data = data.result, 
                 aes(y=RoC, 
                     x= Age))+
      theme_classic()+
      scale_x_continuous(trans = "reverse")+
      coord_flip(ylim=c(0,1))+
      #geom_point(alpha=1/5)+
      geom_line()
    r <- list(ID=dataset.ID, data=data.result, plot=p1)
  }
  
  if (result == "tibble")
  {
    data.result <- data.frame(Age=data.smooth$Age$age[1:sample.size.work], RoC=DC.res.s)
    row.names(data.result) <- row.names(data.smooth$Pollen)[1:sample.size.work]
    r <- data.result
  }
  
  
 end.time <- Sys.time()
 time.length <- end.time - start.time
 print("-")
 print (paste("RATEPOL finished", end.time,"taking",time.length, class(time.length)))
 
 return(r)
 
}

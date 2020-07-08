fc_simulate_pollen_data_in_all_methods <- function(time=0:10e3, 
                                                   nforc=4, 
                                                   mean=100, 
                                                   sdev=.15, 
                                                   nprox=10, 
                                                   var=20, 
                                                   range=15,
                                                   manual.edit = T,
                                                   breaks=c(2000,3000),
                                                   jitter = T,
                                                   rarity=T,
                                                   Working.Unit="levels", 
                                                   BIN.size=500, 
                                                   N.shifts=5, 
                                                   N.datasets=100, 
                                                   interest.treshold=8000){
  
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  performance.list.plot <- vector("list",length = 20);
  
  breaks.seq <- c("empty",rep(c("focus","empty"),length(breaks)))
  
  for(i in 1:N.datasets) {
    
    # create random data
    random.data <- fc_simulate_pollen_data(time=time,
                                           nforc = nforc, 
                                           mean = mean, 
                                           sdev=sdev,
                                           nprox = nprox,
                                           var=var,
                                           range = range,
                                           manual.edit = manual.edit,
                                           breaks = breaks,
                                           jitter = jitter,
                                           rarity = rarity)
    
    for(j in 1:20){
      data.temp<- fc_R_ratepol( data.source.pollen =  random.data$filtered.counts,
                                data.source.age = random.data$list_ages,
                                sm.type = performance.smooth[j],
                                N.points = 5,
                                range.age.max = 500, 
                                grim.N.max = 9,
                                Working.Unit = Working.Unit,
                                BIN.size = BIN.size,
                                N.shifts = N.shifts,
                                rand = 1,
                                standardise = F, 
                                DC = performance.DC[j],
                                interest.treshold = interest.treshold,
                                Peak="GAM",
                                Debug = F) %>%
        as_tibble() %>%
        select(-PEAK)
     
      # PEAK detection 
      
      # Median peak treshold
        # treshold for RoC peaks is set as median of all RoC in dataset
        r.treshold <- median(data.temp$ROC)
        # mark peaks which have 95% quantile above the treshold asPeak.treshold
        data.temp$PEAK.T <- data.temp$ROC.dw>r.treshold
      
      # GAM  
        # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
        pred.gam <-  predict.gam(gam(ROC~s(AGE,k=3), data = data.temp))
        pred.gam.diff <- data.temp$ROC - pred.gam
        data.temp$PEAK.G <- (pred.gam.diff) > 1.5*sd(pred.gam.diff)
      
      # SNI  
        # set moving window of 5 times higher than average distance between samples
        mean.age.window <- 5 * mean( diff(data.temp$AGE) )
        # calculate SNI (singal to noise ratio)
        SNI.calc <- CharSNI(data.frame(data.temp$AGE, data.temp$ROC, pred.gam),mean.age.window)
        # mark points with SNI higher than 3
        data.temp$PEAK.S <- SNI.calc$SNI > 3 & data.temp$ROC > pred.gam
      
        tibble.fin <- data.temp %>%
          mutate(dataset.ID = i,
                 SMOOTH = performance.smooth[j],
                 DC = performance.DC[j]) %>%
          select(dataset.ID, SMOOTH, DC,  Working_Unit ,AGE, ROC, PEAK.T, PEAK.G, PEAK.S)
        
        if (i == 1 & j == 1){
          tibble.fin.sum <- tibble.fin
        } else {
          tibble.fin.sum <- rbind(tibble.fin.sum,tibble.fin)
        }
    }
  }
  return(tibble.fin.sum)
}
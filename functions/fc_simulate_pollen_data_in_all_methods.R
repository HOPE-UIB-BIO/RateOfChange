fc_simulate_pollen_data_in_all_methods <- function(random.data,
                                                   Working.Unit="levels",
                                                   BIN.size=500, 
                                                   N.shifts=5, 
                                                   interest.treshold=8000){
  
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  performance.list.plot <- vector("list",length = 20);
  
  N.datasets <- nrow(random.data)
  
  for(i in 1:N.datasets) {
    
    cat(paste0("dataset ",i,"/",N.datasets), "\n")
    
    for(j in 1:20){
      
      #print (j)
      
      repeat {
        try( data.temp<- fc_estimate_RoC( data_source_pollen =  random.data$filtered.counts[[i]],
                                          data_source_age = random.data$list_ages[[i]],
                                          smooth_method  = performance.smooth[j],
                                          smooth_N_points  = 5,
                                          smooth_age_range  = 500, 
                                          smooth_N_max  = 9,
                                          Working_Units  = Working.Unit,
                                          bin_size  = BIN.size,
                                          Number_of_shifts  = N.shifts,
                                          rand = 1,
                                          treads = 1,
                                          age_uncertainty  = F,
                                          standardise  = F, 
                                          DC  = performance.DC[j],
                                          interest_threshold  = interest.treshold,
                                          Peak_detection ="Threshold",
                                          Debug  = F) %>%
               as_tibble() %>%
               select(-PEAK), silent = T )
        if (exists("data.temp")==T) break
      }
      
     
      # PEAK detection 
      
      # Median peak treshold
        # treshold for RoC peaks is set as median of all RoC in dataset
        r.treshold <- median(data.temp$ROC)
        # mark peaks which have 95% quantile above the treshold asPeak.treshold
        data.temp$PEAK.T <- data.temp$ROC.dw>r.treshold
      
      # GAM  
        # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
        suppressWarnings(try(pred.gam <-   mgcv::predict.gam( mgcv::gam(ROC~s(AGE,k=3), data = data.temp, family = "Gamma",
                                                                        method = "REML"), type="response"),
                             silent =T))
        
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
        
        rm (data.temp)
        
        if (i == 1 & j == 1){
          tibble.fin.sum <- tibble.fin
        } else {
          tibble.fin.sum <- rbind(tibble.fin.sum,tibble.fin)
        }
    }
  }
  return(tibble.fin.sum)
}

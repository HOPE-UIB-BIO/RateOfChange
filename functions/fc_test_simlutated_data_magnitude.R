fc_test_simlutated_data_magnitude <- function(time=0:10e3, 
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
                                             BIN=F, 
                                             BIN.size=500, 
                                             Shiftbin=F, 
                                             N.shifts=5, 
                                             rand.sets=10, 
                                             interest.treshold=8000)
{
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  performance.tibble <- tibble(SMOOTH=NA, DC=NA, RoC_max=NA,RoC_max_SD=NA,RoC_max_05=NA,RoC_max_95=NA, .rows=20);
  
  performance.tibble$SMOOTH <- performance.smooth
  performance.tibble$DC <- performance.DC
  
  for(j in 1:nrow(performance.tibble))
  {
    res.vect <- vector(mode = "double",length = rand.sets)
    
    # for each randomisation
    for(i in 1:rand.sets)
    {
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
                                             rarity = rarity);
      
      data.temp <- fc_ratepol( data.source.pollen =  random.data$filtered.counts,
                              data.source.age = random.data$list_ages,
                              sm.type = performance.smooth[j],
                              N.points = 5,
                              range.age.max = 500, 
                              grim.N.max = 9,
                              BIN = BIN,
                              BIN.size = BIN.size,
                              Shiftbin = Shiftbin,
                              N.shifts = N.shifts,
                              rand = 1,
                              standardise = F, 
                              S.value = 150 ,
                              DC = performance.DC[j],
                              interest.treshold = interest.treshold,
                              Debug = F) %>%
        as.data.frame();
      
      res.vect[i]<- max(data.temp$RUN.RoC)
    }
    
    performance.tibble$RoC_max[j] <- median(res.vect)
    performance.tibble$RoC_max_SD[j] <- sd(res.vect)
    performance.tibble$RoC_max_05[j] <- quantile(res.vect,0.025)
    performance.tibble$RoC_max_95[j] <- quantile(res.vect,0.975)
    
  }
  
  plot.p <- performance.tibble %>%
    ggplot(aes(y=RoC_max, x=SMOOTH, fill=SMOOTH))+
    geom_bar(stat = "identity",position = "dodge",color="gray50")+
    geom_errorbar(aes(ymax=RoC_max_95, ymin=RoC_max_05),position = "dodge",color="gray30", width=0.2, size=0.5)+
    facet_wrap(~DC, scales = "free_y", ncol = 4)+
    theme_classic()+
    theme(strip.background = element_blank(), legend.position = "none")
  
  return(list(data=performance.tibble, plot=plot.p))
}
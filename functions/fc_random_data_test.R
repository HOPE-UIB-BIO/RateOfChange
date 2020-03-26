fc_random_data_test <- function(time=0:10e3, 
                                nforc=4, 
                                mean=100, 
                                sdev=.15, 
                                nprox=10, 
                                var=20, 
                                range=15,
                                manual.edit = T,
                                breaks=c(2000,3000),
                                jitter = T,
                                BIN=F, 
                                BIN.size=500, 
                                Shiftbin=F, 
                                N.shifts=5, 
                                rand.sets=10, 
                                interest.treshold=8000)
{
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  performance.list.plot <- vector("list",length = 20);
    
  breaks.seq <- c("empty",rep(c("focus","empty"),length(breaks)))
    
    for(j in 1:length(performance.list.plot))
    {
      list.res <- data.frame(SEGMENT= character(),SIGNIF=character(), 
                             Value=double(), RAND=double())
      
      # for each randomisation
      for(i in 1:rand.sets)
      {
        # create random data
        random.data <- fc_random_data(time=time,
                                      nforc = nforc, 
                                      mean = mean, 
                                      sdev=sdev,
                                      nprox = nprox,
                                      var=var,
                                      range = range,
                                      manual.edit = manual.edit,
                                      breaks = breaks,
                                      jitter = jitter)
        
        data.temp<- fc_ratepol( data.source.pollen =  random.data$filtered.counts,
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
      
        
        # prelocate space for results
        signif.list <- c("Peak.treshold","Peak.gam","Peak.SNI")
        res.mat <- array(data=NA, dim = c( length(breaks.seq),3), dimnames = list( breaks.seq, signif.list ) )
        DF.res.sig <- data.frame(SEGMENT= character(), DETECT= character(), value=double(),DESC=character(),SIGNIF=character())
        res.mat.temp <- res.mat
        
        #look in arround im the distace of the differences betwen ecollogical breaks
        window.size <- diff(breaks)[1]/2
        
        # extract values for each significance test
        for(l in 1:length(signif.list))
        {
          for (k in 1:length(breaks.seq))
          {
            if(breaks.seq[k]=="empty")
            {
              if(k == 1){
                target <- c(min(time),breaks[k]-window.size)  
              } else{
                if(k==length(breaks.seq))
                {
                  target <- c(breaks[(k-1)/2]+window.size, max(time)) 
                }else{
                  target <- c(breaks[(k-1)/2]+window.size,breaks[(k+1)/2]-window.size)  
                }
              }
            } else {
              target <- c(breaks[k/2]-window.size,breaks[k/2]+window.size)
            }
            
            # sort by values
            target <-sort(target)
            
            # test if there is a peak in selected signif value
            res.mat.temp[k,l] <-  data.temp %>%
              filter(RUN.Age.Pos > target[1] & RUN.Age.Pos < target[2]) %>%
              select(signif.list[l]) %>%
              pull(var=1) %>%
              any()
          }
        }
  
        # save result from single randomization
        sum.temp <-res.mat.temp %>%
          as_tibble() %>%
          mutate(SEGMENT = breaks.seq) %>%
          pivot_longer(-c(SEGMENT)) %>%
          group_by(SEGMENT, name) %>%
          summarise(VALUE = mean(value)) %>% 
          rename(SIGNIF = name)
        
        list.res <- rbind(list.res,data.frame(sum.temp, RAND=i)) 
      }
      
      # summary of randomisation
      plot.data <- list.res %>%
        group_by(SEGMENT,SIGNIF) %>%
        summarise(VALUE.M= median(VALUE, na.rm = T),
                  VALUE.05 = quantile(VALUE,0.025, na.rm = T),
                  VALUE.95 = quantile(VALUE,0.975, na.rm = T)
                  ) %>%
        ungroup()
      
      performance.list.plot[[j]] <- plot.data %>%
        ggplot(aes(y=VALUE.M,x=SEGMENT))+
        geom_bar(stat="identity")+
        geom_errorbar(aes(ymin=VALUE.05, ymax=VALUE.95))+
        facet_wrap(~SIGNIF)+
        theme_classic()+
        #scale_color_manual(values = c("gray30","gray80") )+
        coord_cartesian(ylim = c(0,1))+
        xlab("Detect signal")+ylab("% of detected peaks")+
        ggtitle(paste("smooth",performance.smooth[j],"- DC",performance.DC[j]))
      
    }
    
  p.all <- ggarrange(performance.list.plot[[1]],performance.list.plot[[2]],performance.list.plot[[3]],performance.list.plot[[4]],
                     performance.list.plot[[5]],performance.list.plot[[6]],performance.list.plot[[7]],performance.list.plot[[8]],
                     performance.list.plot[[9]],performance.list.plot[[10]],performance.list.plot[[11]],performance.list.plot[[12]],
                     performance.list.plot[[13]],performance.list.plot[[14]],performance.list.plot[[15]],performance.list.plot[[16]],
                     performance.list.plot[[17]],performance.list.plot[[18]],performance.list.plot[[19]],performance.list.plot[[20]],
                     nrow = 5, ncol = 4, common.legend = T)
  return(p.all)
}
fc_random_data_test <- function(time=0:10e3, 
                                nforc=4, 
                                mean=100, 
                                sdev=.15, 
                                nprox=10, 
                                var=20, 
                                range=15,
                                manual.edit = T,
                                breaks=c(1000,3000,5000,7000),
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
      list.res <- data.frame(SEGMENT= character(), DETECT= character(), 
                             value=double(),DESC=character(),SIGNIF=character(), RAND=double())
      
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
      
        fc_extract_table <- function(x)
        {
          res.temp <- data.temp%>%
            filter(RUN.Age.Pos > target[1] & RUN.Age.Pos < target[2]) %>%
            select(x) %>%
            table()  
          return(res.temp)
        }
        
        # prelocate space for results
        signif.list <- c("Peak.treshold","Peak.gam","Peak.SNI")
        res.mat <- array(data=0, dim = c( length(breaks.seq),2), dimnames = list( breaks.seq, c( FALSE,TRUE) ) )
        DF.res.sig <- data.frame(SEGMENT= character(), DETECT= character(), value=double(),DESC=character(),SIGNIF=character())
        
        #look in arround im the distace of the differences betwen ecollogical breaks
        window.size <- diff(breaks)[1]
        
        # extract values for each significance test
        for(l in 1:length(signif.list))
        {
          res.mat.temp <- res.mat
          
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
            
            target <-sort(target)
            
            res.temp <- fc_extract_table(signif.list[l])
            res.mat.temp[k,1] <- res.temp[1] 
            res.mat.temp[k,2] <- sum(res.temp)-res.temp[1]
          }
          
          sum.temp <-res.mat.temp %>%
            as_tibble() %>%
            mutate(SEGMENT = breaks.seq) %>%
            pivot_longer(-c(SEGMENT)) %>%
            group_by(SEGMENT, name) %>%
            summarise(VALUE = sum (value)) %>%
            pivot_wider(id_cols = c(SEGMENT,name), values_from = VALUE) %>%
            column_to_rownames(var = "SEGMENT")
          
          res<-  c(sum.temp/rowSums(sum.temp))  %>%
             as.data.frame() %>%
             mutate(SEGMENT=c("focus","empty")) %>%
             pivot_longer(-c(SEGMENT)) %>%
             mutate(DESC = c("type II error","right detected signal","right detected space","type I error"),
                    SIGNIF = rep(signif.list[l],4)) %>%
            rename(DETECT = name) 
          
           DF.res.sig<- rbind(DF.res.sig,res)
           
        }
  
        # save result from single randomization
        list.res <- rbind(list.res,data.frame(na.omit(DF.res.sig), RAND=i)) 
      }
      
      # summary of randomisation
      plot.data <- list.res %>%
        group_by(SEGMENT,DETECT,DESC,SIGNIF) %>%
        summarise(VALUE= median(value),
                  VALUE.05 = quantile(value,0.025),
                  VALUE.95 = quantile(value,0.975)
                  ) %>%
        ungroup() %>%
        mutate(RESULT =c(rep("right",3),rep("wrong",6),rep("right",3)))
      
      performance.list.plot[[j]] <- plot.data %>%
        ggplot(aes(y=VALUE,x=SEGMENT, color=RESULT))+
        geom_bar(aes(fill=DESC),stat="identity", position = "dodge")+
        geom_errorbar(aes(ymin=VALUE.05, ymax=VALUE.95, group=DESC),position = "dodge")+
        facet_wrap(~SIGNIF)+
        theme_classic()+
        scale_color_manual(values = c("gray30","gray80") )+
        coord_cartesian(ylim = c(0,1))+
        xlab("Detect signal")+ylab("%")+
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
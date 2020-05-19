fc_test_simlutated_data_succsess <- function(sim.data,
                                             breaks=c(2000,3000)){
  
  performance.smooth <- c(rep("none",4),rep("m.avg",4),rep("grim",4),rep("age.w",4),rep("shep",4));
  performance.DC <- c(rep(c("euc","euc.sd","chord","chisq"),5));
  breaks.seq <- c("empty",rep(c("focus","empty"),length(breaks)))
    
  N.datasets <- sim.data$dataset.ID %>%
    unique() %>%
    length()
  
  pb <- txtProgressBar(min = 0, max = N.datasets, style = 3)
  
  for (i in 1: N.datasets){
    
    setTxtProgressBar(pb, i)
    
    for(j in 1:20){
      
      # prelocate space for results
      signif.list <- c("PEAK.T","PEAK.G","PEAK.S")
      res.mat.temp <- array(data=NA, dim = c( length(breaks.seq),3), dimnames = list( breaks.seq, signif.list ) ) 
      
      #look in arround im the distace of the differences betwen ecollogical breaks
      window.size <- 500 #diff(breaks)[1]/2
      
      # extract values for each PEAK significance test type
      for(l in 1:3){
        
        for (k in 1:length(breaks.seq)){
          
          if(breaks.seq[k]=="empty")
          {
            if(k == 1){
              target <- c(min(sim.data$AGE),breaks[k]-window.size)  
            } else{
              if(k==length(breaks.seq))
              {
                target <- c(breaks[(k-1)/2]+window.size, max(sim.data$AGE)) 
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
          res.mat.temp[k,l] <-  sim.data %>%
            filter(dataset.ID == i) %>%
            filter(SMOOTH ==performance.smooth[j]) %>%
            filter(DC == performance.DC[j]) %>%
            filter(AGE> target[1] & AGE < target[2]) %>%
            select(signif.list[l]) %>%
            pull(var=1) %>%
            any(.,na.rm=T)
        }
        
      }
      
      # save result from single randomization
      sum.temp <- res.mat.temp %>%
        as_tibble() %>%
        mutate(SEGMENT = breaks.seq) %>%
        pivot_longer(-c(SEGMENT)) %>%
        group_by(SEGMENT, name) %>%
        summarise(VALUE = mean(value)) %>%
        mutate(dataset.ID = i,
               SMOOTH =performance.smooth[j],
               DC = performance.DC[j]) %>%
        rename(PEAK = name) %>%
        dplyr::select(dataset.ID, SMOOTH, DC, PEAK, SEGMENT, VALUE)
      
      if (i ==1 & j ==1 ){
        sum.temp.res <- sum.temp 
      } else {
        sum.temp.res <- rbind(sum.temp.res,sum.temp)
      }
    }
  }
  
  close(pb)
  
  # summary of randomisation
  SUM.data <- sum.temp.res %>%
    group_by(SMOOTH,DC,PEAK,SEGMENT) %>%
    summarise(VALUE.M = mean(VALUE, na.rm = T),
              VALUE.SD = sd(VALUE, na.rm = T),
              VALUE.05 = quantile(VALUE,0.025, na.rm = T),
              VALUE.95 = quantile(VALUE,0.975, na.rm = T)
    ) %>%
    ungroup()
  
  return(list(RawData =sum.temp.res, SumData =SUM.data))    
}

fc_ratepol <- function (data.source, standardise = T, S.value = 150, sm.type = "age.w", N.points = 5, range.age.max = 300, grim.N.max = 9)
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
  
  # data extraction (already include data check)
  data.work <- fc_extract(data.sub)
  
  # standardisation of pollen data to X(S.value) number of pollen grains 
  if(standardise==T) # 
  {
    data.work$Pollen <- fc_standar(data.work$Pollen, S.value)
    
    if(any(rowSums(data.work$Pollen)!=S.value))
      stop("standardisation was unsuccesfull")
  }
  
  # data check with proportioning
  data.work <- fc_check(data.work, proportion = T)
  
  # smoothing of pollen data
  data.smooth <- fc_smooth(data.work, 
                           sm.type = sm.type, 
                           N.points = N.points,
                           grim.N.max = grim.N.max, 
                           range.age.max = range.age.max)
  #data check with proportioning
  data.smooth <- fc_check(data.smooth, proportion = T)
  
  
  # calculate DC for each sample
  # exploration
  data.source <- data.smooth
  rm(data.source)
  test1<- fc_calDC(data.smooth,DC="euc")
  test1
  test2<- fc_calDC(data.smooth,DC="euc.sd")
  test2
  
  cor(test1, test2)
  
  p1<-ggplot(data = reshape2::melt(data.smooth$Pollen), 
             aes(y=value, 
                 x=c(rep(1:nrow(data.smooth$Pollen),ncol(data.smooth$Pollen)) )))+
    theme_classic()+
    coord_flip(ylim = c(0,1))+
    geom_point(alpha=1/5)+
    geom_line(group=1)+
    facet_wrap(~variable)
  
  p2<-ggplot(data=data.frame(value=test,age=1:length(test)), aes(y=age, x=value))+geom_point()+theme_classic()
  
  plot(p1)
  plot(p2)
  
  # rest of the function
  # WIP
  
}
fc_draw_RoC <- function (data.source, type="map", age.treshold = 15000, dataset.N="" )
{
  
  res.df.plot <- data.source %>%
    select(dataset.id, collection.handle, long, lat, ROC) %>%
    unnest(cols = c(ROC)) %>%
    select(.,-c(newage))
  
  if (type=="perplot")
  {
    p.fin<- res.df.plot %>%
      ggplot(aes( y=RoC.median, 
                  x= age))+
      theme_classic()+
      scale_x_continuous(trans = "reverse")+
      coord_flip(xlim=c(0,age.treshold), ylim = c(0,5))+
      geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), alpha=1/5)+
      geom_line(alpha=1, size=1)+
      geom_point(data = res.df.plot[res.df.plot$Peak==T,],color="blue", alpha=1, size=1)+
      geom_hline(yintercept = 0, color="red")+
      facet_wrap(~ dataset.id)+
      xlab("Age")+ylab("Rate of Change")
  }
  
  if( type=="summary")
  {
    RoC_summary_p1 <- res.df.plot %>%
      filter (age < age.treshold) %>%
      ggplot(aes( y=RoC.median, 
                  x= age))+
      theme_classic()+
      scale_x_continuous(trans = "reverse")+
      coord_flip(xlim=c(0,age.treshold), ylim = c(0,5))+
      geom_line(aes(group=as.factor(dataset.id)),alpha=1/10, size=1)+
      geom_hline(yintercept = 0, color="red")+
      geom_smooth(color="green", method = "loess", se=F)+
      xlab("Age")+ylab("Rate of Change")

    RoC_summary_p1b <-res.df.plot %>%
      filter (age < age.treshold) %>%
      ggplot(aes(x=age))+
      geom_density(fill="gray")+
      theme_classic()+
      scale_x_continuous(trans = "reverse")+
      coord_flip(xlim=c(0,age.treshold))+
      xlab("")+ylab("Density of the samples")
    
    RoC_summary_p2 <- res.df.plot %>%
      filter (age < age.treshold) %>%
      filter(Peak==T) %>%
      ggplot(aes( y=RoC.median, 
                  x= age))+
      theme_classic()+
      scale_color_gradient2(low="white",mid="darkblue",high="black", midpoint = 4)+
      scale_x_continuous(trans = "reverse")+
      coord_flip(xlim=c(0,age.treshold), ylim = c(0,5))+
      geom_point(aes(color=RoC.median), alpha=1/10, size=3)+
      geom_hline(yintercept = 0, color="red")+
      geom_smooth(color="orange", method = "loess", se=F)+
      xlab("")+ylab("Rate of Change")+
      theme(legend.position = "none")

    RoC_summary_p2b <- res.df.plot %>%
      filter (age < age.treshold) %>%
      filter (Peak==T) %>%
      ggplot(aes( x= age))+
      theme_classic()+
      scale_x_continuous(trans = "reverse")+
      coord_flip(xlim=c(0,age.treshold))+
      geom_density(fill="gray")+
      xlab("")+ylab("Density of Peak-points")
    
    
    
    p.fin <- ggarrange(RoC_summary_p1,RoC_summary_p1b,
                            RoC_summary_p2,RoC_summary_p2b,
                            ncol = 4, labels = c("A","B","C","D"))
  }
  
  if (type=="map")
  {
    lat.dim <- c(min(data.source$lat),max(data.source$lat)) 
    long.dim <- c(min(data.source$long),max(data.source$long))
    
    p.fin <- data.source %>%
      mutate(
        N.RoC.points = select(.,ROC) %>%
          pluck(.,1) %>% 
          map_dbl(.,.f=function(x) {
            dplyr::filter(x,Peak==T) %>%
              nrow() }),
        N.Samples = select(.,ROC) %>%
          pluck(.,1) %>% 
          map_dbl(.,.f=function(x) {
            nrow(x) }),
        Ratio = N.RoC.points/N.Samples
      ) %>% 
      ggplot(aes(x = long, y = lat)) +
      borders(fill = "gray60", colour = "gray50") +
      coord_fixed(ylim = lat.dim, xlim = long.dim) +
      geom_point(aes(color=Ratio, size=N.RoC.points)) + 
      scale_color_gradient("Ration of Peak-points to total samples",low="black",high = "red")+
      scale_size( "Number of Peak-points")+
      labs(x = "Longitude", y = "Latitude")+
      theme_classic()
  }
  
  if(type=="singleplot")
  {
    if(is.numeric(dataset.N)==F)
      stop("Dataset number is not a valid number")
    
    single.plot <- data.source %>%
      filter(dataset.id==dataset.N) %>%
      select(ROC) %>%
      unnest(cols = ROC)
    
    p.fin <- single.plot %>%
      ggplot(aes( y=RoC.median, 
                  x= age)) +
      theme_classic() +
      scale_x_continuous(trans = "reverse") +
      coord_flip(xlim=c(0,age.treshold), ylim = c(0,5)) +
      geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), alpha=1/5) +
      geom_line(alpha=1, size=2) +
      geom_point(data = . %>% filter(Peak==T),color="blue", alpha=1, size=3) +
      geom_hline(yintercept = single.plot$RoC.median %>%
                   median(), color="blue") +
      geom_hline(yintercept = 0, color="red") +
      xlab("Age")+ylab("Rate of Change")
  }
  
 return(p.fin) 
}
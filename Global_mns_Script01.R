test <- tibble_Europe_Roc %>%
  select(dataset.id, collection.handle, long, lat, ROC) %>%
  unnest(cols = c(ROC)) %>%
  mutate(BIN = ceiling(RUN.Age.Pos/500))
  

DF.sum <- data.frame(matrix(nrow = 16, ncol=7))
names(DF.sum) <- c("BIN","N.sample","N.peaks","RoC","RoC.05","RoC.95","Peak.R")

for(i in 1:length(unique(test$BIN)))
{
  samples = 1000
  
    DF.sum$BIN[i] = i
  
    X<- test %>%
      filter(BIN == i)  
    
    DF.sum$N.sample[i] = nrow(X)
    
    X.peak <- test %>%
      filter(BIN == i & Peak.gam == T)
    
    
    DF.sum$N.peaks[i] = nrow(X.peak)
    
    X.sample <- sample_n(X,size = samples, replace = T)
    
    DF.sum$RoC[i] = median(X.sample$RUN.RoC)
    DF.sum$RoC.05[i] <- quantile(X.sample$RUN.RoC, 0.025)
    DF.sum$RoC.95[i] <- quantile(X.sample$RUN.RoC, 0.975)

    DF.sum$Peak.R[i] = ((X.sample %>%
                       filter(Peak.gam == T) %>%
                       nrow())/samples) *100 
        
}


BIN.names <- c("0-500","500-1000","1500-2000","2000-2500","2500-3000","3000-3500","3500-4000",
               "4000-4500","4500-5000","5000-5500","5500-6000","6000-6500","6500-7000",
               "7000-7500","7500-8000","8000-8500")

ggarrange(
  test %>%
    ggplot(aes(x= BIN))+
    geom_density(fill="gray80",color="gray50")+
    theme_classic()+
    scale_x_continuous(name= "time",
                       breaks = DF.sum$BIN,
                       labels = BIN.names)+
    ylab("density of samples"),
DF.sum %>%
  ggplot(aes(y=RoC, x= BIN))+
  geom_jitter(data=test, aes(x=BIN, y=RUN.RoC),color="gray50", alpha=1/10)+
  geom_point(aes(y=Peak.R/20, x= BIN), color="blue")+
  geom_line(aes(y=Peak.R/20, x= BIN), color="blue")+
  geom_point(aes(size=Peak.R), shape=15)+
  geom_errorbar(aes(ymin=RoC.05,ymax=RoC.95), width=0.2)+
  geom_line()+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_continuous(name= "time",
                     breaks = DF.sum$BIN,
                     labels = BIN.names)+
  scale_y_continuous(sec.axis = sec_axis(~ .*20, name = "% of Peak points"))+
  coord_cartesian(ylim=c(0,1.5)),
nrow = 2, align = "v"
)



ggsave("Test.res.plot.pdf",
       width = 35,height = 20, units = "cm")


ggarrange(
test %>%
  ggplot(aes(y =  lat, x= RUN.Age.Pos))+
  geom_point(aes(size=RUN.RoC, color=RUN.RoC), alpha=1/5)+
  theme_classic()+
  coord_cartesian()+
  scale_color_gradient2(low = "white",mid = "deepskyblue", high = "darkblue", midpoint = 1)+
  xlab("time")+ylab("latitude"),
test %>%
  group_by(lat) %>%
  summarise(LAT=mean(lat)) %>%
  ggplot(aes(x=LAT))+
  geom_density(fill="gray80",color="gray50")+
  theme_classic()+
  coord_flip()+
  xlab("latitude"),
ncol = 2, align = "h")


ggsave("Test.res.plot2.pdf",
       width = 35,height = 20, units = "cm")





test$lat.BIN<- floor(test$lat/10)*10  

lat.BIN.list<- test$lat.BIN %>%
  unique %>%
  sort


DF.lat.sum <- data.frame(matrix(nrow = 80, ncol = 4))
names(DF.lat.sum) <- c("lat.BIN","BIN","RoC","Var")

for(i in 1:length(unique(test$lat.BIN)))
{
  for(k in 1:length(unique(test$BIN)))
  {
    print(paste(i,k))
    DF.lat.sum$lat.BIN[((i-1)*16) + k] = lat.BIN.list[i]
    DF.lat.sum$BIN[((i-1)*16) + k] = k
    
    X<- test %>%
      filter(lat.BIN == lat.BIN.list[i] & BIN == k)  
    
    X.sample <- sample_n(X,size = 1000, replace = T)
  
    DF.lat.sum$RoC[((i-1)*16) + k] = median(X.sample$RUN.RoC)
    DF.lat.sum$Var[((i-1)*16) + k] = var(X.sample$RUN.RoC)
      
  }
}


DF.lat.sum %>%
  ggplot(aes(y=lat.BIN,x=BIN))+
  geom_point(aes(size=RoC, color=Var))+
  scale_color_gradient2(low = "deepskyblue", mid = "darkblue", high = "black", midpoint = 0.04)+
  scale_size(range = c(1,50),guide = "none")+
  theme_classic()+
  scale_x_continuous(name= "time",
                     breaks = DF.sum$BIN,
                     labels = BIN.names)+
  ylab("latitude")

ggsave("Test.res.plot3.pdf",
       width = 35,height = 20, units = "cm")

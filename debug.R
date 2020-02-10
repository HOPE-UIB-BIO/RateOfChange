data.source <- tibble_Europe2[2,]
rand = 99
standardise = T
S.value = 150
sm.type = "grim" 
N.points = 5
range.age.max = 300
grim.N.max = 9
DC = "chisq"
Debug = F

test <- fc_ratepol(data.sub[2,],
                  rand = 99,
                  standardise = T, 
                  S.value = 150, 
                  sm.type = "m.avg", 
                  N.points = 5, 
                  range.age.max = 300, 
                  grim.N.max = 9,
                  DC = "chord",
                  Debug = F)

test$Data %>% ggplot(aes( y=RoC.mean, 
                     x= DF.Age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), color="gray")+
  #geom_point(alpha=1/5)+
  geom_line()+
  geom_point(data = test$Data[test$Data$Peak==T,], color="red", size=3)+
  geom_hline(yintercept = median(test$Data$RoC.mean), color="blue")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip()

data.sub<-tibble_Europe2[c(1:10),]


data.frame(POLLEN=reshape2::melt(tibble_Europe2[2,]$filtered.counts[[1]]), 
           AGE=rep(tibble_Europe2[2,]$list_ages[[1]]$ages$age, ncol(tibble_Europe2[2,]$filtered.counts[[1]]))) %>%
  ggplot(aes( y=POLLEN.value, 
                  x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  #geom_point(alpha=1/5)+
  #geom_line(aes(group=POLLEN.variable), alpha=1/10)+
  #geom_smooth(method = "loess",color="blue",se=F)+
  geom_density_2d()+
  xlab("Age")+ylab("Pollen")+
  coord_flip()



n.sets <- 1:length(test)

names(n.sets) <- n.sets

for(i in n.sets[c(1:65,67:70,72:74,76:125,127:131,133:length(n.sets))])
{
  vec<- test[[i]]$age
  for(j in 2:length(vec))
    {
    z<- abs(vec[j]-vec[j-1])
    if(z==0) stop(names(n.sets)[i])
    }
}

# 66, 71, 75, 126, 132

tibble_Europe2$list_ages[[66]]$ages
tibble_Europe2$list_ages[[75]]$ages
tibble_Europe2$list_ages[[126]]$ages

tibble_Europe2$list_ages[[71]]$ages
tibble_Europe2$list_ages[[132]]$ages




r.m %>% ggplot(aes( y=RoC.mean, 
                     x= DF.Age))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RoC.05q, ymax=RoC.95q), color="gray")+
  #geom_point(alpha=1/5)+
  geom_line()+
  geom_point(data = r.m[r.m$Peak==T,], color="red", size=3)+
  geom_hline(yintercept = median(r.m$RoC.mean), color="blue")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip()
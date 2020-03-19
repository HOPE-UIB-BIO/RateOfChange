which(tibble_Europe2$dataset.id %in%  25318 )

dataset.25318

dataset.N <- 68
data.source.pollen <- tibble_Europe2$filtered.counts[[dataset.N]]
data.source.age <- tibble_Europe2$list_ages[[dataset.N]]
sm.type = "age.w" 
N.points = 5
range.age.max = 500 
grim.N.max = 9
BIN = T
BIN.size = 500
Shiftbin  = T
N.shifts = 5
rand = 100
standardise = T 
S.value = 150 
DC = "chisq"
interest.treshold = 8000
Debug = F

data.source.extrap <- data.bin


data.source <- data.sd

data.subset<- data.work
BINS <- SELECTED.BINS

data.source.bin<- data.smooth.check

data.source.check <- list.res

data.source <- data.smooth


data.source.pollen.extract <- data.source.pollen
data.source.age.extract <- data.source.age

data.source.pollen.check <- data.source.pollen
data.source.age.check <- data.source.age

ages = r.m.full$RUN.Age.Pos;
CHAR = r.m.full$RUN.RoC;
thresh = pred.gam;
BandWidth = mean.age.window


ggarrange(
  data.frame(ROC=r.m.full$RUN.RoC, AGE=r.m.full$RUN.Age.Pos, GAM=pred.gam) %>%
    ggplot()+
    geom_line(aes(y=ROC, x=AGE))+
    geom_line(aes(y=GAM,x=AGE))+
    geom_point(aes(y=ROC,x=AGE))+
    geom_point(data=r.m.full[r.m.full$Peak.SNI == T,], aes(y=RUN.RoC,x=RUN.Age.Pos), color="red"),
  data.frame(SNI=SNI.calc$SNI, AGE=r.m.full$RUN.Age.Pos) %>%
    ggplot()+
    geom_line(aes(y=SNI, x=AGE, group=1))+
    geom_point(aes(y=SNI, x=AGE))+
    geom_hline(yintercept = 3),
  nrow = 2
)


dataset.N <- 2

test <- fc_ratepol( data.source.pollen =  tibble_Europe2$filtered.counts[[dataset.N]],
                    data.source.age = tibble_Europe2$list_ages[[dataset.N]],
                    sm.type = "age.w" ,
                    N.points = 5,
                    range.age.max = 500, 
                    grim.N.max = 9,
                    BIN = T,
                    BIN.size = 500,
                    Shiftbin = T,
                    N.shifts = 5,
                    rand = 1000,
                    standardise = T, 
                    S.value = 150 ,
                    DC = "chisq",
                    interest.treshold = 8000,
                    Debug = F)


test %>% ggplot(aes( y=RUN.RoC, 
                     x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
  geom_line(size=1)+
  geom_line(data=data.frame(RUN.RoC = predict.gam(gam(RUN.RoC~s(RUN.Age.Pos), data = test)),
                            RUN.Age.Pos = test$RUN.Age.Pos),
            color="blue", size=1)+
  geom_point(color="black",size=1)+
  geom_point(data = test[test$Peak.treshold==T,], color="yellow", size=1)+
  geom_point(data = test[test$Peak.treshold.95==T,], color="orange", size=2)+
  geom_point(data = test[test$Peak.gam==T,], color="red", size=3)+
  geom_point(data = test[test$Peak.SNI==T,], color="purple", size=4)+
  geom_hline(yintercept = median(test$RUN.RoC), color="green")+
  geom_hline(yintercept = 0, color="red")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,5))





test %>% ggplot(aes( y=RUN.RoC, 
                     x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
  #geom_point(alpha=1/5)+
  geom_line(size=1)+
  #geom_point(data = test[test$soft.Peak==T,], color="orange", size=3)+
  #geom_point(data = test[test$Peak==T,], color="blue", size=3)+
  #geom_hline(yintercept = mean(test$RoC.median), color="green")+
  #geom_hline(yintercept = median(test$RUN.RoC), color="blue")+
  geom_hline(yintercept = 0, color="red")+
  geom_smooth(method = "loess")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,5))

example <- fc_extract(tibble_Europe2$filtered.counts[[dataset.N]],
           tibble_Europe2$list_ages[[dataset.N]],
           Debug=F) %>%
  fc_smooth(., 
            sm.type = "shep", 
            N.points = 5,
            range.age.max = 500,
            grim.N.max = 9,
            Debug=F) %>%
  fc_check(.,proportion = T)


N.taxa <- 10

Common.list <- example$Pollen %>%
  colSums() %>% 
  sort(decreasing = T) %>%
  .subset(.,1:N.taxa) %>%
  names()


data.frame(POLLEN=example$Pollen %>%
             select(.,Common.list) %>%
             reshape2::melt() , 
           AGE=rep(example$Age$newage, N.taxa))  %>%
  ggplot(aes( y=POLLEN.value, 
              x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=rep(0,length(POLLEN.value)),ymax=POLLEN.value, fill=POLLEN.variable), 
              color="gray20", alpha=1/5, size=0.1)+
  #geom_smooth(method = "loess",color="blue",se=F)+
  #facet_wrap(~POLLEN.variable, ncol = N.taxa)+
  xlab("Age")+ylab("Pollen")+
  coord_flip(xlim=c(0,8000), ylim = c(0,1))

ggsave("smooth_shep.png")



Common.list <- tibble_Europe2$filtered.counts[[dataset.N]] %>%
  colSums() %>% 
  sort(decreasing = T) %>%
  .subset(.,1:N.taxa) %>%
  names()


data.frame(POLLEN=tibble_Europe2$filtered.counts[[dataset.N]] %>%
                 select(.,Common.list) %>%
                 reshape2::melt() , 
               AGE=rep(tibble_Europe2$list_ages[[dataset.N]]$ages$age, N.taxa)) %>%
  ggplot(aes( y=POLLEN.value, 
              x= AGE))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=rep(0,length(POLLEN.value)),ymax=POLLEN.value), color="gray20", fill="gray80")+
  facet_wrap(~POLLEN.variable)+
  xlab("Age")+ylab("Pollen")+
  coord_flip(xlim=c(0,8000), ylim = c(0,1000))



s.time <- Sys.time()

tibble_Europe_Roc.test <-  dataset.25318 %>%
  mutate(., ROC = map2(filtered.counts,list_ages,
                       .f = function(.x,.y)
                       {res <- fc_ratepol(
                         data.source.pollen = .x,
                         data.source.age = .y,
                         sm.type = "age.w", 
                         N.points = 5, 
                         range.age.max = 500, 
                         grim.N.max = 9,
                         BIN = T,
                         BIN.size = 500,
                         Shiftbin = T,
                         N.shifts = 5,
                         rand = 1000,
                         standardise = T, 
                         S.value = 150, 
                         DC = "chisq",
                         interest.treshold = 8000,
                         Debug = F
                       )} ))

f.time <- Sys.time()
tot.time <- f.time - s.time
tot.time


fc_draw_RoC(tibble_Europe_Roc.test,type = "singleplot", 
            dataset.N = 25318, age.treshold = 8000)

r.m.full %>% ggplot(aes( y=RUN.RoC, 
                      x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
  #geom_point(alpha=1/5)+
  geom_line(size=1)+
  #geom_point(data = r.m.full[r.m.full$soft.Peak==T,], color="blue", size=3)+
  #geom_hline(yintercept = mean(test$RoC.median), color="green")+
  geom_hline(yintercept = median(r.m.full$RUN.RoC), color="blue")+
  geom_hline(yintercept = 0, color="red")+
  xlab("Age")+ylab("Rate of Change")+
  coord_flip(xlim=c(0,8000), ylim=c(0,3)) 


r.m.full %>% ggplot(aes( x=RUN.Age.Pos))+
  geom_density(fill="gray80")+
  theme_classic()+
  coord_flip()+
  scale_x_continuous(trans = "reverse")
  
  


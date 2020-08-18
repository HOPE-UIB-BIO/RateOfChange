which(tibble_Europe$dataset.id %in%  17334 )

dataset.25318

dataset.N <- 2
data.source.pollen <- tibble_Europe$filtered.counts[[dataset.N]]
data.source.age <- tibble_Europe$list_ages[[dataset.N]]
sm.type = "shep" 
N.points = 5
range.age.max = 500 
grim.N.max = 9
Working.Unit = "MW"
BIN.size = 500
N.shifts = 5
rand = 99
standardise = T 
S.value = 150 
DC = "chord"
Peak = "GAM"
interest.treshold = 8000
Debug = F



time=seq(from=0, to=10e3, by=100);
time= tibble_Europe2$list_ages[[2]]$ages$age
nforc=4;
mean=100; 
sdev=.15; 
nprox=10; 
var=20;
range=15;
manual.edit = T;
breaks=c(2000,3000);
jitter = T;
rarity=T;
BIN.size=500; 
N.shifts=5;
Working.Unit="MW";
N.datasets=100;
interest.treshold=8000;


random.data = sim_ld_recent


data.source.pollen =  random.data$filtered.counts;
data.source.age = random.data$list_ages;
sm.type = performance.smooth[j];
N.points = 5;
range.age.max = 500; 
grim.N.max = 9;
Working.Unit = Working.Unit;
BIN.size = BIN.size;
N.shifts = N.shifts;
rand = 1;
standardise = F; 
DC = performance.DC[j];
interest.treshold = interest.treshold;
Peak="GAM";
Debug = F;












rm(test)
test <- fc_random_data_test(time= tibble_Europe2$list_ages[[2]]$ages$age,
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

ggsave("Perfomance_test_exp.pdf", units = "cm",width = 30, height = 20)

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


data.temp %>% ggplot(aes( y=RUN.RoC, 
                     x= RUN.Age.Pos))+
  theme_classic()+
  scale_x_continuous(trans = "reverse")+
  geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
  geom_line(size=1)+
  geom_line(data=data.frame(RUN.RoC = predict.gam(gam(RUN.RoC~s(RUN.Age.Pos), data = data.temp)),
                            RUN.Age.Pos = data.temp$RUN.Age.Pos),
            color="blue", size=1)+
  geom_point(color="black",size=1)+
  geom_point(data = data.temp[data.temp$Peak.treshold==T,], color="yellow", size=1)+
  geom_point(data = data.temp[data.temp$Peak.treshold.95==T,], color="orange", size=2)+
  geom_point(data = data.temp[data.temp$Peak.gam==T,], color="red", size=3)+
  geom_point(data = data.temp[data.temp$Peak.SNI==T,], color="purple", size=4)+
  geom_hline(yintercept = median(data.temp$RUN.RoC), color="green")+
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
  
  
random.data <- fc_random_data(time=tibble_Europe2$list_ages[[2]]$ages$age,
                              nforc=4,
                              mean=100, 
                              sdev=.15, 
                              nprox=10, 
                              var=20,
                              range=15,
                              manual.edit = T,
                              breaks=c(2000,2500,5000,5500),
                              jitter = T
                              )


random.data.ex <- fc_extract(random.data$filtered.counts, random.data$list_ages) %>%
  fc_smooth("none") %>%
  fc_check(., proportion = T)

Roc.random <- fc_ratepol(data.source.pollen= random.data$filtered.counts, 
           data.source.age= random.data$list_ages,
           sm.type = "none",
           standardise = F,
           BIN=F,
           BIN.size=500, 
           Shiftbin=F,
           N.shifts=5,
           rand = 1,
           interest.treshold=8000,
           DC="euc"
           )

breaks=c(2000,2500,5000,5500);
ggarrange(as.data.frame(random.data.ex$Pollen) %>%
            mutate(AGE = random.data.ex$Age$newage) %>%
            pivot_longer(., cols = -c(AGE)) %>%
            ggplot(aes(x=AGE, y=value))+
            geom_ribbon(aes(ymin=rep(0,length(value)), ymax=value, fill=name), alpha=1/5, color="gray30")+
            geom_vline(xintercept = breaks, color="red")+
            coord_flip(ylim = c(0,1),xlim=c(0,8000))+
            facet_wrap(~name, ncol=nprox)+
            theme_classic()+
            scale_x_continuous(trans = "reverse")+
            theme(legend.position ="none"),
          Roc.random %>% ggplot(aes( y=RUN.RoC, 
                               x= RUN.Age.Pos))+
            theme_classic()+
            scale_x_continuous(trans = "reverse")+
            geom_ribbon(aes(ymin=RUN.RoC.05q, ymax=RUN.RoC.95q), color="gray", alpha=1/5)+
            geom_line(size=1)+
            geom_line(data=data.frame(RUN.RoC = predict.gam(gam(RUN.RoC~s(RUN.Age.Pos), data = Roc.random)),
                                      RUN.Age.Pos = Roc.random$RUN.Age.Pos),
                      color="blue", size=1)+
            geom_point(color="black",size=1)+
            geom_point(data = Roc.random[Roc.random$Peak.treshold==T,], color="yellow", size=1)+
            geom_point(data = Roc.random[Roc.random$Peak.treshold.95==T,], color="orange", size=2)+
            geom_point(data = Roc.random[Roc.random$Peak.gam==T,], color="red", size=3)+
            geom_point(data = Roc.random[Roc.random$Peak.SNI==T,], color="purple", size=4)+
            geom_hline(yintercept = median(Roc.random$RUN.RoC), color="green")+
            geom_hline(yintercept = 0, color="red")+
            geom_vline(xintercept = breaks, color="red")+
            geom_vline(xintercept = c(breaks-500,breaks+500), color="gray80", lty=2)+
            xlab("Age")+ylab("Rate of Change")+
            coord_flip(xlim=c(0,8000), ylim=c(0,2)),
          ncol=2
)




test<- fc_ratepol(data.source.pollen = data.source.pollen,
                  data.source.age = data.source.age,
                  standardise = F,
                  rand=1,
                  BIN = F,
                  Shiftbin = F,
                  sm.type = "none",
                  DC="euc",
                  interest.treshold = 8000)


ggarrange(as.data.frame(forcing) %>%
            mutate(AGE = time) %>%
            pivot_longer(., cols = -c(AGE)) %>%
            ggplot(aes(y=AGE, x= value))+
            geom_point(aes(color=name))+
            geom_hline(yintercept = breaks, color="red")+
            theme_classic()+
            coord_cartesian(ylim=c(8000,0))+
            scale_y_continuous(trans = "reverse")+
            theme(legend.position = "none"),
          fc_extract(data.source.pollen, data.source.age) %>%
            fc_smooth("none") %>%
            fc_check(., proportion = T) %>%
            pluck("Pollen") %>%
            mutate(AGE = time) %>%
            pivot_longer(-c(AGE)) %>%
            ggplot(aes(x=AGE, y=value))+
            geom_ribbon(aes(ymin=rep(0,length(value)), ymax=value, fill=name), alpha=1/5, color="gray30")+
            geom_vline(xintercept = breaks, color="red")+
            coord_flip(xlim=c(0,8000), ylim=c(0,1))+
            #facet_wrap(~name, ncol=nprox)+
            theme_classic()+
            scale_x_continuous(trans = "reverse")+
            ylab("% of pollen grains")+
            theme(legend.position = "none"),
          test %>%
            ggplot(aes( y=RUN.RoC, 
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
            geom_vline(xintercept = breaks, color="red")+
            xlab("Age")+ylab("Rate of Change")+
            coord_flip(xlim=c(0,8000), ylim=c(0,1)),
          ncol = 3)

ggsave("RandomData_exp2.pdf", units = "cm", width = 30, height = 20)
          



  
